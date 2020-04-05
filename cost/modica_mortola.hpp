//
// Created by janos on 27.11.19.
//

#pragma once


#include <ceres/first_order_function.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <memory>

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

