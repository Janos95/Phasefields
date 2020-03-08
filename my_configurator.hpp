//
// Created by janos on 07.12.19.
//


#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <unsupported/Eigen/CXX11/Tensor>

struct ConfigurationTraits
{
    using RealType = double;
    using ScalarType = Eigen::Matrix<RealType, 1, 1>;
    using VectorType = Eigen::VectorXd;
    using SparseMatrixType = Eigen::SparseMatrix<double, Eigen::RowMajor, int>;
    using FullMatrixType = Eigen::MatrixXd;

    using VecType = Eigen::Matrix<RealType, 3, 1>;
    using MatType = Eigen::Matrix<RealType, 3, 3>;

    using TripletType = Eigen::Triplet<double>;
    using MaskType = Eigen::Matrix<bool, Eigen::Dynamic, 1>;

    using TensorType = Eigen::Tensor< RealType, 3>;
};
