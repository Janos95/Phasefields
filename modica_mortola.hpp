//
// Created by janos on 20.11.19.
//

#pragma once

#include "cost/quadrature_ref_triangle.hpp"

#include <igl/grad.h>
#include <igl/doublearea.h>
#include <igl/massmatrix.h>
#include <igl/massmatrix_intrinsic.h>

#include <Eigen/Sparse>
#include <Eigen/Geometry>

#include <vector>



auto interpolateBasisFunction(const double w1, const double w2, const int basis)
{
    double interp;
    switch(basis)
    {
        case 0: interp = (1. - w1 - w2);
            break;
        case 1: interp = w1;
            break;
        case 2: interp = w2;
            break;
        default:
            assert(false);
            __builtin_unreachable();
    }
    return interp;
}



struct DoubleWellPotential
{
        double v1, v2, v3;

        auto operator()(const double w1, const double w2) const
        {
            auto vInterp =  w1*v2 + w2 * v3 + (1. - w1 - w2) * v1;
            return 9./16. * std::pow(vInterp * vInterp - 1., 2);
        }
};

struct DoubleWellPotentialGrad
{
    double v1, v2, v3;
    int i;

    auto operator()(const double w1, const double w2) const
    {
        auto vInterp =  w1*v2 + w2 * v3 + (1. - w1 - w2) * v1;
        auto pInterp = interpolateBasisFunction(w1, w2, i);

        return 9./4. * ( vInterp * vInterp - 1.) * vInterp * pInterp;
    }
};


struct IndicatorFunction
{
    double v1, v2, v3;

    auto operator()(const double w1, const double w2) const
    {
        double vInterp =  w1*v2 + w2 * v3 + (1. - w1 - w2) * v1;
        return 1./4. * std::pow( vInterp + 1., 2);
    }
};

struct IndicatorFunctionGrad
{
    double v1, v2, v3;
    int i;

    auto operator()(const double w1, const double w2) const
    {
        auto vInterp = w1*v2 + w2 * v3 + (1. - w1 - w2) * v1;
        auto pInterp = interpolateBasisFunction(w1, w2, i);

        return .5 * ( vInterp - 1.) * pInterp;
    }
};


class ModicaMortola
{
public:

    ModicaMortola(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const double epsilon):
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

    double apply (const Eigen::VectorXd &U) const
    {
        return .5 * (
                m_epsilon * U.transpose() * m_GSQ * U +
                doubleWellPotential(U) / m_epsilon);
    }

    Eigen::VectorXd derivative(const Eigen::VectorXd &U) const
    {
        //TODO here there are some cashing opportunities...
        return m_epsilon * m_GSQ * U + doubleWellGradient(U);
    }

private:

    double doubleWellPotential(const Eigen::VectorXd& U) const
    {
        QuadratureRefTriangle<double> quad;
        double value(0);

        for (int i = 0; i < m_F.cols(); ++i) {
            auto f = m_F.col(i);
            DoubleWellPotential pot{U[f[0]],U[f[1]],U[f[2]]};
            value += m_dblA[i] * quad.integrate(pot);
        }

        return value;
    }

    Eigen::VectorXd doubleWellGradient(const Eigen::VectorXd& U) const
    {
        QuadratureRefTriangle<double> quad;
        Eigen::VectorXd dwG = Eigen::VectorXd::Zero(U.rows());

        for (int i = 0; i < m_F.cols(); ++i) {
            auto f = m_F.col(i);
            DoubleWellPotentialGrad pot{U[f[0]],U[f[1]],U[f[2]]};

            for (int j = 0; j < 3; ++j) {
                pot.i = j;
                dwG[f[j]] += m_dblA[i] * quad.integrate(pot);
            }
        }

        return dwG;
    }

    const Eigen::MatrixXd& m_V;
    const Eigen::MatrixXi& m_F;

    Eigen::SparseMatrix<double> m_GSQ;
    Eigen::VectorXd m_dblA;

    double m_epsilon;
};


class PhaseArea
{

public:

    PhaseArea(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F):
            m_F(F),
            m_V(V),
            m_dblA(F.rows())
    {
        igl::doublearea(m_V, m_F, m_dblA);
    }

    double apply(const Eigen::VectorXd& U) const
    {
        QuadratureRefTriangle<double> quad;
        double value = 0;

        for (int i = 0; i < m_F.cols(); ++i) {
            auto f = m_F.col(i);
            IndicatorFunction ind{U[f[0]],U[f[1]],U[f[2]]};
            value += m_dblA[i] * quad.integrate(ind);
        }

        return value;
    }

    Eigen::VectorXd derivative(const Eigen::VectorXd &U) const
    {
        QuadratureRefTriangle<double> quad;
        Eigen::VectorXd iG = Eigen::VectorXd::Zero(U.rows());


        for (int i = 0; i < m_F.cols(); ++i) {
            auto f = m_F.col(i);
            IndicatorFunctionGrad interp{U[f[0]],U[f[1]],U[f[2]]};
            for (int j = 0; j < 3; ++j) {
                interp.i = j;
                iG[f[j]] += m_dblA[i] * quad.integrate(interp);
            }
        }

        return iG;
    }

private:

    const Eigen::MatrixXd& m_V;
    const Eigen::MatrixXi& m_F;

    Eigen::VectorXd m_dblA;
};
