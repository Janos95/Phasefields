//
// Created by janos on 05.03.20.
//

#include <connectedness_constraint.hpp>
#include <phasefield_initialization.hpp>

#include <Magnum/Primitives/Capsule.h>
#include <Magnum/Primitives/Grid.h>
#include <Magnum/EigenIntegration/Integration.h>
#include <Magnum/Trade/MeshData.h>
#include <Magnum/Math/Vector2.h>
#include <Magnum/Math/Vector3.h>
#include <Magnum/Magnum.h>

#include <Corrade/TestSuite/Tester.h>

using namespace Corrade;
using namespace Magnum;
#include <iostream>
namespace {


    struct ConnectednessConstraintTest: TestSuite::Tester {
        Eigen::VectorXd m_U;
        Eigen::MatrixXi m_F;
        Eigen::MatrixXd m_V;
        Eigen::VectorXd m_grad;
        ConnectednessConstraint m_auto,m_ana;
        double m_r;
        
        ConnectednessConstraintTest()
        {
            std::tie(m_V, m_F) = toEigen(Primitives::capsule3DSolid(20, 20, 20, 1.));
            m_U.resize(m_V.rows());
            initRandomNormal(m_U);
            m_U *= 100;
            m_auto = ConnectednessConstraint(m_V, m_F, 1e-2, .85, .95, GradientFlag::Automatic);
            m_ana = ConnectednessConstraint(m_V, m_F, 1e-2, .85, .95, GradientFlag::Analytic);
            m_grad = Eigen::VectorXd(m_V.rows());
            addTests({&ConnectednessConstraintTest::grid, &ConnectednessConstraintTest::capsule});
            addBenchmarks({&ConnectednessConstraintTest::benchAutodiff, &ConnectednessConstraintTest::benchAnalytic}, 10);
        }

        void testMesh(Trade::MeshData const& mesh){
            auto [V, F] = toEigen(mesh);

            Eigen::VectorXd U(V.rows());
            initRandomNormal(U);
            U *= 100;
            double epsilon = 1e-2;

            ConnectednessConstraint constraint(V, F, epsilon, .85, .95, GradientFlag::Analytic);
            ConnectednessConstraint constraintAutodiff(V, F, epsilon, .85, .95, GradientFlag::Automatic);

            double r1,r2;
            Eigen::VectorXd grad1(V.rows());
            Eigen::VectorXd grad2(V.rows());
            constraint.Evaluate(U.data(), &r1, grad1.data());
            constraintAutodiff.Evaluate(U.data(), &r2, grad2.data());

            auto n = std::count_if(grad1.begin(), grad1.end(), [](const auto& x){ return std::abs(x) < std::numeric_limits<double>::epsilon(); });
            auto m = std::count_if(grad2.begin(), grad2.end(), [](const auto& x){ return std::abs(x) < std::numeric_limits<double>::epsilon(); });

            CORRADE_VERIFY(n==m);
            CORRADE_VERIFY(std::abs(r2 - r1) < 1e-8);
            CORRADE_VERIFY((grad1 - grad2).norm() < 1e-6);
        }

        void grid(){
            testMesh(Primitives::grid3DSolid({20,20}));
            testMesh(Primitives::grid3DSolid({10,10}));
            testMesh(Primitives::grid3DSolid({30,30}));
        }

        void capsule(){
            testMesh(Primitives::capsule3DSolid(10, 10, 10, 1));
            testMesh(Primitives::capsule3DSolid(5, 5, 5, 1));
        }

        void benchAutodiff(){
            m_auto.Evaluate(m_U.data(), &m_r, m_grad.data());
        }
        
        void benchAnalytic(){
            m_ana.Evaluate(m_U.data(), &m_r, m_grad.data());
        }
        
        };

}

CORRADE_TEST_MAIN(ConnectednessConstraintTest)
