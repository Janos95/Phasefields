//
// Created by janos on 27.11.19.
//

#pragma once

#include "interpolation.hpp"
#include "quadrature_ref_triangle.hpp"

#include <igl/grad.h>
#include <igl/doublearea.h>

#include <ceres/dynamic_cost_function.h>
#include <ceres/first_order_function.h>
#include <ceres/loss_function.h>

#include <Corrade/Utility/Assert.h>
#include <Corrade/TestSuite/Comparator.h>

#include <memory>



class SumProblem : public ceres::FirstOrderFunction, public ceres::CostFunction {
public:
    using ProblemVecType = std::vector<std::unique_ptr<ceres::FirstOrderFunction>>;
    using LossFunctionVecType = std::vector<std::unique_ptr<ceres::LossFunction>>;
    using MappedVectorType = Eigen::Map<Eigen::VectorXd>;


    explicit SumProblem(int problemSize);


    void push_back(
            std::unique_ptr<ceres::FirstOrderFunction> func,
            std::unique_ptr<ceres::LossFunction> loss = nullptr,
            double weight = 1);

    std::vector<double> computeSeperateCosts(const Eigen::VectorXd& U);

    bool Evaluate(const double* parameters,
                  double* cost,
                  double* jacobians) const override;


    bool Evaluate(double const* const* parameters,
                          double* residuals,
                          double** jacobians) const override;


    int NumParameters() const override;

private:
    ProblemVecType m_problems;
    LossFunctionVecType m_losses;
    LossFunctionVecType m_owner;
    std::vector<double> m_weights;
};

class InterfaceEnergy : public ceres::FirstOrderFunction
{
public:

    InterfaceEnergy(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const double epsilon);

    int NumParameters() const override;

    bool Evaluate(double const* parameters,
                  double* residual,
                  double* jacobians) const override;

private:

    double m_epsilon;
    const Eigen::MatrixXd& m_V;
    const Eigen::MatrixXi& m_F;

    Eigen::SparseMatrix<double> m_GSQ;
    Eigen::VectorXd m_dblA;

};


class PotentialEnergy : public ceres::FirstOrderFunction
{
public:

    PotentialEnergy(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const double epsilon);

    int NumParameters() const override;

    bool Evaluate(double const* parameters,
                  double* cost,
                  double* jacobians) const override;

private:

    double m_epsilon;
    const Eigen::MatrixXd& m_V;
    const Eigen::MatrixXi& m_F;

    Eigen::VectorXd m_dblA;

};


class AreaRegularizer : public ceres::FirstOrderFunction
{

public:

    AreaRegularizer(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

    int NumParameters() const override;

    bool Evaluate(double const* parameters,
                  double* cost,
                  double* jacobians) const override;

private:

    const Eigen::MatrixXd& m_V;
    const Eigen::MatrixXi& m_F;

    Eigen::VectorXd m_dblA;
    double m_area;
};

