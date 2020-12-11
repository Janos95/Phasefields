//
// Created by janos on 12/5/20.
//

#include "CudaSolver.h"
#include "CuSparse/cusparse_cholesky_solver.h"

namespace Phasefield {

CuSparseCholeskySolver<double>* cast(void* p) {
    return (CuSparseCholeskySolver<double>*)p;
}

void CUDASolver::compute(Eigen::SparseMatrix<double>& A) {
    auto solverPtr = CuSparseCholeskySolver<double>::create(A.cols());
    m_solver = solverPtr.release();
    m_matrix = A;

    if (m_doOrdering) {
        // compute permutation
        //PermutationMatrix P;
        //Ordering ordering;
        //ordering(Acsr.selfadjointView<Eigen::Upper>(), P);

        // set permutation to solver
        //solver->setPermutaion(n, P.indices().data());
    }

    auto solver = cast(m_solver);
    solver->analyze(m_matrix.nonZeros(), m_matrix.outerIndexPtr(), m_matrix.innerIndexPtr());
    solver->factorize(m_matrix.valuePtr());
}

Eigen::VectorXd CUDASolver::solve(Eigen::VectorXd const& b) {
    Eigen::VectorXd x(b.size());
    cast(m_solver)->solve(b.data(), x.data());
    return x;
}

Eigen::ComputationInfo CUDASolver::info() {
    switch(cast(m_solver)->info()) {
        case CuSparseCholeskySolver<double>::SUCCESS : return Eigen::Success;
        default : return Eigen::NumericalIssue;
    }
}

CUDASolver::~CUDASolver() {
    delete (CuSparseCholeskySolver<double>*)m_solver;
}

}