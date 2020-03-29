//
// Created by janos on 2/23/20.
//
#include "modica_mortola.hpp"
#include "quadrature_ref_triangle.hpp"
#include "interpolation.hpp"

#include <Corrade/Utility/Assert.h>

#include <igl/doublearea.h>
#include <igl/grad.h>

InterfaceEnergy::InterfaceEnergy(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const double epsilon):
        m_epsilon(epsilon),
        m_V(V),
        m_F(F),
        m_dblA(F.rows())
{
    igl::doublearea(m_V, m_F, m_dblA);
    CORRADE_INTERNAL_ASSERT(!m_dblA.hasNaN());

    Eigen::SparseMatrix<double> G;
    igl::grad(m_V,m_F,G);
    G.makeCompressed();

#ifndef NODEBUG
    Eigen::Map<Eigen::VectorXd> map(G.valuePtr(), G.nonZeros());
    CORRADE_INTERNAL_ASSERT(!map.hasNaN());
#endif

    const auto & T = (m_dblA.replicate(3,1) * .5).asDiagonal();
    m_GSQ = G.transpose() * T * G;

}


int InterfaceEnergy::NumParameters() const {
    return m_V.rows();
}


bool InterfaceEnergy::Evaluate(double const* parameters,
              double* residual,
              double* jacobians) const
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

PotentialEnergy::PotentialEnergy(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const double epsilon):
        m_epsilon(epsilon),
        m_V(V),
        m_F(F),
        m_dblA(F.rows())
{
    igl::doublearea(m_V, m_F, m_dblA);
}

int PotentialEnergy::NumParameters() const {
    return m_V.rows();
}

bool PotentialEnergy::Evaluate(double const* parameters,
              double* cost,
              double* jacobians) const
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



AreaRegularizer::AreaRegularizer(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F):
        m_V(V),
        m_F(F),
        m_dblA(F.rows())
{
    igl::doublearea(m_V, m_F, m_dblA);
    m_area = m_dblA.sum() / 2.;
}

int AreaRegularizer::NumParameters() const {
    return m_V.rows();
}


bool AreaRegularizer::Evaluate(double const* parameters,
              double* cost,
              double* jacobians) const
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


