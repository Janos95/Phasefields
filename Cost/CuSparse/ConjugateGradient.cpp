//
// Created by janos on 12/7/20.
//

#include "ConjugateGradient.h"

namespace Phasefield {

void ConjugateGradient::solve() {
    /* Get handle to the CUBLAS context */
    cublasHandle_t cublasHandle = 0;
    cublasStatus_t cublasStatus;
    cublasStatus = cublasCreate(&cublasHandle);

    checkCudaErrors(cublasStatus);

    /* Get handle to the CUSPARSE context */
    cusparseHandle_t cusparseHandle = 0;
    checkCudaErrors(cusparseCreate(&cusparseHandle));

    checkCudaErrors(cudaMalloc((void **)&d_col, nz*sizeof(int)));
    checkCudaErrors(cudaMalloc((void **)&d_row, (N+1)*sizeof(int)));
    checkCudaErrors(cudaMalloc((void **)&d_val, nz*sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&d_x, N*sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&d_r, N*sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&d_p, N*sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&d_Ax, N*sizeof(float)));

    /* Wrap raw data into cuSPARSE generic API objects */
    cusparseSpMatDescr_t matA = NULL;
    checkCudaErrors(cusparseCreateCsr(
            &matA, N, N, nz, d_row, d_col, d_val, CUSPARSE_INDEX_32I,
            CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F));
    cusparseDnVecDescr_t vecx = NULL;
    checkCudaErrors(cusparseCreateDnVec(&vecx, N, d_x, CUDA_R_32F));
    cusparseDnVecDescr_t vecp = NULL;
    checkCudaErrors(cusparseCreateDnVec(&vecp, N, d_p, CUDA_R_32F));
    cusparseDnVecDescr_t vecAx = NULL;
    checkCudaErrors(cusparseCreateDnVec(&vecAx, N, d_Ax, CUDA_R_32F));

    /* Initialize problem data */
    cudaMemcpy(d_col, J, nz*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_row, I, (N+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_val, val, nz*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_x, x, N*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_r, rhs, N*sizeof(float), cudaMemcpyHostToDevice);

    alpha = 1.0;
    alpham1 = -1.0;
    beta = 0.0;
    r0 = 0.;

    /* Allocate workspace for cuSPARSE */
    size_t bufferSize = 0;
    checkCudaErrors(cusparseSpMV_bufferSize(
            cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecx,
            &beta, vecAx, CUDA_R_32F, CUSPARSE_MV_ALG_DEFAULT, &bufferSize));
    void *buffer = NULL;
    checkCudaErrors(cudaMalloc(&buffer, bufferSize));

    /* Begin CG */
    checkCudaErrors(cusparseSpMV(
            cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecx,
            &beta, vecAx, CUDA_R_32F, CUSPARSE_MV_ALG_DEFAULT, &buffer));

    cublasSaxpy(cublasHandle, N, &alpham1, d_Ax, 1, d_r, 1);
    cublasStatus = cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);

    k = 1;

    while (r1 > tol*tol && k <= max_iter)
    {
        if (k > 1)
        {
            b = r1 / r0;
            cublasStatus = cublasSscal(cublasHandle, N, &b, d_p, 1);
            cublasStatus = cublasSaxpy(cublasHandle, N, &alpha, d_r, 1, d_p, 1);
        }
        else
        {
            cublasStatus = cublasScopy(cublasHandle, N, d_r, 1, d_p, 1);
        }

        checkCudaErrors(cusparseSpMV(
                cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA,
                vecp, &beta, vecAx, CUDA_R_32F, CUSPARSE_MV_ALG_DEFAULT, &buffer));
        cublasStatus = cublasSdot(cublasHandle, N, d_p, 1, d_Ax, 1, &dot);
        a = r1 / dot;

        cublasStatus = cublasSaxpy(cublasHandle, N, &a, d_p, 1, d_x, 1);
        na = -a;
        cublasStatus = cublasSaxpy(cublasHandle, N, &na, d_Ax, 1, d_r, 1);

        r0 = r1;
        cublasStatus = cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
        cudaDeviceSynchronize();
        printf("iteration = %3d, residual = %e\n", k, sqrt(r1));
        k++;
    }

    cudaMemcpy(x, d_x, N*sizeof(float), cudaMemcpyDeviceToHost);

    float rsum, diff, err = 0.0;

    for (int i = 0; i < N; i++)
    {
        rsum = 0.0;

        for (int j = I[i]; j < I[i+1]; j++)
        {
            rsum += val[j]*x[J[j]];
        }

        diff = fabs(rsum - rhs[i]);

        if (diff > err)
        {
            err = diff;
        }
    }

    cusparseDestroy(cusparseHandle);
    cublasDestroy(cublasHandle);
    if (matA       ) { checkCudaErrors(cusparseDestroySpMat(matA)); }
    if (vecx       ) { checkCudaErrors(cusparseDestroyDnVec(vecx)); }
    if (vecAx      ) { checkCudaErrors(cusparseDestroyDnVec(vecAx)); }
    if (vecp       ) { checkCudaErrors(cusparseDestroyDnVec(vecp)); }

    free(I);
    free(J);
    free(val);
    free(x);
    free(rhs);
    cudaFree(d_col);
    cudaFree(d_row);
    cudaFree(d_val);
    cudaFree(d_x);
    cudaFree(d_r);
    cudaFree(d_p);
    cudaFree(d_Ax);
}


}