//
// Created by janos on 12/14/20.
//

#include "SparseMatrix.h"


SparseMatrix::SparseMatrix(size_t m, size_t n, size_t nnz):
        L(*this)
{
    data = cholmod_l_spzeros(m, n, nnz, CHOLMOD_REAL, common);
}

 SparseMatrix::SparseMatrix(Triplets& T):
        L(*this)
{
    cholmod_triplet *triplet = T.toCholmod();
    data = cholmod_l_triplet_to_sparse(triplet, triplet->nnz, common);
}

SparseMatrix::~SparseMatrix()
{
    cholmod_l_free_sparse(&data, common);
    data = NULL;
}

SparseMatrix SparseMatrix::identity(size_t m, size_t n)
{
    return SparseMatrix(cholmod_l_speye(m, n, CHOLMOD_REAL, common));
}

SparseMatrix SparseMatrix::transpose() const
{
    return SparseMatrix(cholmod_l_transpose(data, 1, common));
}

size_t SparseMatrix::nRows() const
{
    return data->nrow;
}

size_t SparseMatrix::nCols() const
{
    return data->ncol;
}

size_t SparseMatrix::nnz() const
{
    return cholmod_l_nnz(data, common);
}

double SparseMatrix::norm(int norm) const
{
    return cholmod_l_norm_sparse(data, norm, common);
}

cholmod_sparse* SparseMatrix::toCholmod()
{
    return data;
}

 void scale(double s, cholmod_sparse *A)
{
    // A = s*A
    DenseMatrix S(1, 1);
    S(0, 0) = s;
    cholmod_l_scale(S.toCholmod(), CHOLMOD_SCALAR, A, common);
}

 cholmod_sparse* add(cholmod_sparse *A, cholmod_sparse *B, double alpha[2], double beta[2])
{
    // C = alpha*A + beta*B
    return cholmod_l_add(A, B, alpha, beta, 1, 1, common);
}

 cholmod_sparse* mul(cholmod_sparse *A, cholmod_sparse *B)
{
    // C = A*B
    return cholmod_l_ssmult(A, B, 0, 1, 1, common);
}

 void mul(cholmod_sparse *A, cholmod_dense *X, cholmod_dense *Y, double alpha[2], double beta[2])
{
    // Y = alpha*(A*X) + beta*Y
    cholmod_l_sdmult(A, 0, alpha, beta, X, Y, common);
}

 SparseMatrix operator*(const SparseMatrix& A, double s)
{
    cholmod_sparse *data = A.copy();
    scale(s, data);

    return SparseMatrix(data);
}

 SparseMatrix operator+(const SparseMatrix& A, const SparseMatrix& B)
{
    double alpha[2] = {1.0, 1.0};
    double beta[2] = {1.0, 1.0};
    return SparseMatrix(add(A.data, B.data, alpha, beta));
}

 SparseMatrix operator-(const SparseMatrix& A, const SparseMatrix& B)
{
    double alpha[2] = {1.0, 1.0};
    double beta[2] = {-1.0, -1.0};
    return SparseMatrix(add(A.data, B.data, alpha, beta));
}

 SparseMatrix operator*(const SparseMatrix& A, const SparseMatrix& B)
{
    return SparseMatrix(mul(A.data, B.data));
}

 DenseMatrix operator*(const SparseMatrix& A, const DenseMatrix& X)
{
    DenseMatrix Y(A.nRows(), X.nCols());
    double alpha[2] = {1.0, 1.0};
    double beta[2] = {0.0, 0.0};
    mul(A.data, X.data, Y.data, alpha, beta);

    return Y;
}

 SparseMatrix& operator*=(SparseMatrix& A, double s)
{
    scale(s, A.data);
    A.L.clearNumeric();

    return A;
}

 SparseMatrix& operator+=(SparseMatrix& A, const SparseMatrix& B)
{
    double alpha[2] = {1.0, 1.0};
    double beta[2] = {1.0, 1.0};
    A = add(A.data, B.data, alpha, beta);

    return A;
}

 SparseMatrix& operator-=(SparseMatrix& A, const SparseMatrix& B)
{
    double alpha[2] = {1.0, 1.0};
    double beta[2] = {-1.0, -1.0};
    A = add(A.data, B.data, alpha, beta);

    return A;
}

 void SparseMatrix::clear()
{

}

 Triplets::Triplets(size_t m_, size_t n_):
        m(m_),
        n(n_),
        capacity(m_)
{
    data = cholmod_l_allocate_triplet(m, n, capacity, 0, CHOLMOD_REAL, common);
    data->nnz = 0;
}

 Triplets::~Triplets()
{
    clear();
}

 void Triplets::add(size_t i, size_t j, double x)
{
    if (data->nnz == capacity) increaseCapacity();

    ((size_t *)data->i)[data->nnz] = i;
    ((size_t *)data->j)[data->nnz] = j;
    ((double *)data->x)[data->nnz] = x;
    data->nnz++;
}

cholmod_triplet* Triplets::toCholmod()
{
    return data;
}

void Triplets::increaseCapacity()
{
    // create triplet with increased capacity
    capacity *= 2;
    cholmod_triplet *newData = cholmod_l_allocate_triplet(m, n, capacity, 0, CHOLMOD_REAL, common);
    memcpy(newData->i, data->i, data->nzmax*sizeof(size_t));
    memcpy(newData->j, data->j, data->nzmax*sizeof(size_t));
    memcpy(newData->x, data->x, data->nzmax*sizeof(double));
    newData->nnz = data->nnz;

    // clear old triplet and assign the newly created triplet to it
    clear();
    data = newData;
}

void Triplets::clear()
{
    cholmod_l_free_triplet(&data, common);
    data = NULL;
}
