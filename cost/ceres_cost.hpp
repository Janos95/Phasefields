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


    explicit SumProblem(int problemSize)
    {
        mutable_parameter_block_sizes()->push_back(problemSize);
        set_num_residuals(1);
    }


    void push_back(
            std::unique_ptr<ceres::FirstOrderFunction> func,
            std::unique_ptr<ceres::LossFunction> loss = nullptr,
            double weight = 1){
        m_problems.push_back(std::move(func));
        if(!loss || std::abs(weight - 1) < 1e-6)
            m_losses.push_back(std::move(loss));
        else{
            m_losses.push_back(std::make_unique<ceres::ScaledLoss>(loss.get(), weight, ceres::Ownership::DO_NOT_TAKE_OWNERSHIP));
            m_owner.push_back(std::move(loss));
        }
    }

    std::vector<double> computeSeperateCosts(const Eigen::VectorXd& U){
        std::vector<double> costs;
        for(const auto& pb: m_problems){
            double cost = 0;
            pb->Evaluate(U.data(), &cost, nullptr);
            costs.push_back(cost);
        }
        return costs;
    }

    bool Evaluate(const double* parameters,
                  double* cost,
                  double* jacobians) const override {
        auto n = NumParameters();
        *cost = 0.;
        if(jacobians)
            std::fill_n(jacobians, n, 0.);

        auto singleJac = jacobians ? new double[n] : nullptr;

        for(std::size_t i = 0; i < m_problems.size(); ++i){
            double residual = 0;
            m_problems[i]->Evaluate(parameters, &residual, singleJac);
            double out[3] = {residual, 1.};
            if(m_losses[i])
                m_losses[i]->Evaluate(residual, out);
            if(jacobians){
                MappedVectorType mappedJac(jacobians, n), mappedSingleJac(singleJac, n);
                mappedJac.noalias() += out[1] * mappedSingleJac;
            }
            *cost += out[0];
        }

        delete[] singleJac;

        return true;
    }


    bool Evaluate(double const* const* parameters,
                          double* residuals,
                          double** jacobians) const override
    {
        auto jacs = jacobians && *jacobians ? *jacobians : nullptr;
        return Evaluate(*parameters, residuals, jacs);
    }


    int NumParameters() const override{
        auto numParams = m_problems.front()->NumParameters();
        for(const auto& prob: m_problems)
            assert(prob->NumParameters() == numParams);
        return numParams;
    }
private:
    ProblemVecType m_problems;
    LossFunctionVecType m_losses;
    LossFunctionVecType m_owner;
    std::vector<double> m_weights;
};

class InterfaceEnergy : public ceres::FirstOrderFunction
{
public:

    InterfaceEnergy(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const double epsilon):
            m_epsilon(epsilon),
            m_F(F),
            m_V(V),
            m_dblA(F.rows())
    {
        igl::doublearea(m_V, m_F, m_dblA);

        Eigen::SparseMatrix<double> G;
        igl::grad(m_V,m_F,G);
        const auto & T = (m_dblA.replicate(3,1) * .5).asDiagonal();
        m_GSQ = G.transpose() * T * G;
    }


    int NumParameters() const override{
        return m_V.rows();
    }


    bool Evaluate(double const* parameters,
                  double* residual,
                  double* jacobians) const override
    {
        Eigen::Map<const Eigen::VectorXd> U(parameters, m_V.rows()); //TODO: dont const cast

        Eigen::VectorXd intermResult = m_GSQ * U;
        residual[0] = .5 * m_epsilon * U.transpose() * intermResult;

        if(jacobians)
        {
            for (int i = 0; i < intermResult.rows(); ++i) {
                jacobians[i] = m_epsilon * intermResult[i];
            }
        }

        return true;
    }

private:

    const Eigen::MatrixXd& m_V;
    const Eigen::MatrixXi& m_F;

    Eigen::SparseMatrix<double> m_GSQ;
    Eigen::VectorXd m_dblA;

    double m_epsilon;
};


class PotentialEnergy : public ceres::FirstOrderFunction
{
public:

    PotentialEnergy(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const double epsilon):
            m_epsilon(epsilon),
            m_F(F),
            m_V(V),
            m_dblA(F.rows())
    {
        igl::doublearea(m_V, m_F, m_dblA);
    }

    int NumParameters() const override{
        return m_V.rows();
    }

    bool Evaluate(double const* parameters,
                  double* cost,
                  double* jacobians) const override
    {
        Eigen::Map<const Eigen::VectorXd> U(parameters, m_V.rows()); //TODO: dont const cast

        if(jacobians)
            std::fill(jacobians, jacobians + U.rows(), 0);

        QuadratureRefTriangle<double> quad;
        double residual = 0;

        for (int i = 0; i < m_F.rows(); ++i) {
            auto f = m_F.row(i);
            DoubleWellPotential<double> pot{U[f[0]],U[f[1]],U[f[2]]};
            residual += m_dblA[i] * quad.integrate(pot);

            if(jacobians)
            {
                DoubleWellPotentialGrad<double> potGrad{U[f[0]],U[f[1]],U[f[2]]};
                for (int j = 0; j < 3; ++j) {
                    potGrad.i = j;
                    jacobians[f[j]] += m_dblA[i] * quad.integrate(potGrad);
                }
            }
        }

        cost[0] = .5 * residual / m_epsilon;

        if(jacobians)
        {
            Eigen::Map<Eigen::VectorXd> J(jacobians, m_V.rows());
            J *= 1./(2. * m_epsilon);
        }

        return true;
    }

private:

    const Eigen::MatrixXd& m_V;
    const Eigen::MatrixXi& m_F;

    Eigen::VectorXd m_dblA;

    double m_epsilon;
};


class AreaRegularizer : public ceres::FirstOrderFunction
{

public:

    AreaRegularizer(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F):
            m_F(F),
            m_V(V),
            m_dblA(F.rows())
    {
        igl::doublearea(m_V, m_F, m_dblA);
        m_area = m_dblA.sum() / 2.;
    }

    int NumParameters() const override{
        return m_V.rows();
    }


    bool Evaluate(double const* parameters,
                  double* cost,
                  double* jacobians) const override
    {
        Eigen::Map<const Eigen::VectorXd> U(parameters, m_V.rows()); //TODO: dont const cast

        if(jacobians)
            std::fill(jacobians, jacobians + U.rows(), 0);

        QuadratureRefTriangle<double> quad;
        double residual = 0;

        for (int i = 0; i < m_F.rows(); ++i) {
            auto f = m_F.row(i);
            IndicatorFunction<double> ind{U[f[0]],U[f[1]],U[f[2]]};
            residual += m_dblA[i] * quad.integrate(ind);

            if(jacobians)
            {
                IndicatorFunctionGrad<double> interp{U[f[0]],U[f[1]],U[f[2]]};
                for (int j = 0; j < 3; ++j) {
                    interp.i = j;
                    jacobians[f[j]] += m_dblA[i] * quad.integrate(interp);
                }
            }
        }

        //printf("AreaRegularizer: phasefield area, half total area=(%f,%f)\n", residual, m_area/2);
        cost[0] = std::pow(residual - m_area / 2., 2) / 2.;

        if(jacobians)
            Eigen::Map<Eigen::VectorXd>(jacobians, m_V.rows()) *= residual - m_area / 2.;

        return true;
    }


private:

    const Eigen::MatrixXd& m_V;
    const Eigen::MatrixXi& m_F;

    Eigen::VectorXd m_dblA;
    double m_area;
};

