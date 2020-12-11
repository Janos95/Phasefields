//
// Created by janos on 12/5/20.
//

#pragma once

#include <Eigen/SparseCore>

namespace Phasefield {

class CUDASolver {
public:
    explicit CUDASolver() = default;
    ~CUDASolver();

    CUDASolver(const CUDASolver&) = delete;
    CUDASolver& operator=(const CUDASolver&) = delete;

    void compute(Eigen::SparseMatrix<double>&);
    Eigen::ComputationInfo info();

    Eigen::VectorXd solve(Eigen::VectorXd const&);
private:
    Eigen::SparseMatrix<double, Eigen::StorageOptions::RowMajor> m_matrix;
    bool m_doOrdering = false;
    void* m_solver = nullptr;
};

}


