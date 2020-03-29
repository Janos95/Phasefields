//
// Created by janos on 28.11.19.
//

#include <phasefield_initialization.hpp>
#include <modica_mortola.hpp>


#include <ceres/gradient_checker.h>
#include <ceres/dynamic_cost_function.h>

#include <Eigen/Core>

#include <Magnum/Primitives/Grid.h>
#include <Corrade/TestSuite/Tester.h>
#include <Magnum/Trade/MeshData.h>
#include <Magnum/EigenIntegration/Integration.h>
#include <Magnum/Math/Vector2.h>
#include <Magnum/Math/Vector3.h>
#include <Magnum/Magnum.h>

using namespace Corrade;
using namespace Magnum;
namespace {


    constexpr double epsilon = 1e-3;
    struct ModicaMortolaTest : TestSuite::Tester {

        ModicaMortolaTest()
        {
            auto grid = Primitives::grid3DSolid({5,5});

            auto varr = grid.positions3DAsArray();
            auto farr = grid.indicesAsArray();

            using MatrixXU = Eigen::Matrix<UnsignedInt, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
            F = Eigen::Map<const MatrixXU>(farr.data(), farr.size() / 3, 3).cast<int>();
            V = Eigen::MatrixXd (varr.size(), 3);
            for (std::size_t i = 0; i < grid.vertexCount(); ++i) {
                V.row(i) = EigenIntegration::cast<Eigen::Vector3f>(varr[i]).cast<double>();
            }
            U.resize(V.rows());
            initRandomNormal(U);
            U *= 100;

            addTests({
                &ModicaMortolaTest::simple,
                &ModicaMortolaTest::scaled,
                &ModicaMortolaTest::scaledAndWeighted,
                &ModicaMortolaTest::withArea});
        }


        void simple(){
            auto interfaceLoss = std::make_unique<ceres::CauchyLoss>(1);
            auto potentialLoss = std::make_unique<ceres::CauchyLoss>(1);
            auto* cost = new SumProblem(U.size());
            cost->push_back(std::make_unique<InterfaceEnergy>(V,F,epsilon), std::move(interfaceLoss));
            cost->push_back(std::make_unique<PotentialEnergy>(V,F,epsilon), std::move(potentialLoss));
            ceres::GradientChecker gradient_checker(cost, nullptr, ceres::NumericDiffOptions{});
            ceres::GradientChecker::ProbeResults results;
            auto data = U.data();
            CORRADE_VERIFY(gradient_checker.Probe(&data, 1e-3, &results));

        }

        void scaled(){
            auto interfaceLoss = std::make_unique<ceres::CauchyLoss>(11);
            auto potentialLoss = std::make_unique<ceres::CauchyLoss>(34);
            auto* cost = new SumProblem(U.size());
            cost->push_back(std::make_unique<InterfaceEnergy>(V,F,epsilon), std::move(interfaceLoss));
            cost->push_back(std::make_unique<PotentialEnergy>(V,F,epsilon), std::move(potentialLoss));
            ceres::GradientChecker gradient_checker(cost, nullptr, ceres::NumericDiffOptions{});
            ceres::GradientChecker::ProbeResults results;
            auto data = U.data();
            CORRADE_VERIFY(gradient_checker.Probe(&data, 1e-3, &results));
        }

        void scaledAndWeighted(){
            auto interfaceLoss = std::make_unique<ceres::CauchyLoss>(63);
            auto potentialLoss = std::make_unique<ceres::CauchyLoss>(100);
            auto* cost = new SumProblem(U.size());
            cost->push_back(std::make_unique<InterfaceEnergy>(V,F,epsilon), std::move(interfaceLoss), 42.5);
            cost->push_back(std::make_unique<PotentialEnergy>(V,F,epsilon), std::move(potentialLoss), 33.1);
            ceres::GradientChecker gradient_checker(cost, nullptr, ceres::NumericDiffOptions{});
            ceres::GradientChecker::ProbeResults results;
            auto data = U.data();
            CORRADE_VERIFY(gradient_checker.Probe(&data, 1e-3, &results));
        }


        void withArea(){
            auto interfaceLoss = std::make_unique<ceres::CauchyLoss>(1);
            auto potentialLoss = std::make_unique<ceres::CauchyLoss>(2);
            auto areaLoss = std::make_unique<ceres::CauchyLoss>(3);
            auto* cost = new SumProblem(U.size());
            cost->push_back(std::make_unique<InterfaceEnergy>(V,F,epsilon), std::move(interfaceLoss));
            cost->push_back(std::make_unique<PotentialEnergy>(V,F,epsilon), std::move(potentialLoss));
            cost->push_back(std::make_unique<AreaRegularizer>(V,F), std::move(areaLoss), 123);
            ceres::GradientChecker gradient_checker(cost, nullptr, ceres::NumericDiffOptions{});
            ceres::GradientChecker::ProbeResults results;
            auto data = U.data();
            CORRADE_VERIFY(gradient_checker.Probe(&data, 1e-3, &results));
        }

        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        Eigen::VectorXd U;
    };

}

CORRADE_TEST_MAIN(ModicaMortolaTest)
