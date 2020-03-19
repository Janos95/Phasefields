

#include "ceres_cost.hpp"
#include "color_callback.hpp"
#include "phasefield_initialization.hpp"
#include "connectedness_constraint.hpp"
#include "load_mesh_data.hpp"
#include "viewer.hpp"
#include "upload.hpp"
#include "write_mesh.hpp"
#include "colormaps.hpp"
#include "draw_callbacks.hpp"

#include <scoped_timer/scoped_timer.hpp>

#include <Magnum/Trade/MeshData.h>
#include <Magnum/Primitives/Plane.h>
#include <Magnum/Primitives/Cube.h>
#include <Magnum/Primitives/Capsule.h>
#include <Magnum/Primitives/Grid.h>
#include <Magnum/Primitives/Icosphere.h>
#include <Magnum/Image.h>
#include <Magnum/PixelFormat.h>
#include <Magnum/MeshTools/Interleave.h>
#include <Magnum/MeshTools/Transform.h>
#include <Magnum/GL/TextureFormat.h>
#include <Magnum/ImageView.h>

#include <Corrade/Utility/Configuration.h>

#include <ceres/problem.h>
#include <ceres/solver.h>
#include <ceres/loss_function.h>
#include <ceres/gradient_problem.h>
#include <ceres/gradient_problem_solver.h>

#include <fmt/core.h>

#include <thread>
#include <algorithm>
#include <iostream>

using namespace Magnum;
using namespace Corrade;
using namespace Corrade::Containers;

using namespace std::chrono_literals;

auto viewerTest(Trade::MeshData& meshdata)
{
    int argc = 0;
    Viewer viewer(argc, nullptr);
    auto& scene = viewer.scene;
    //viewer.scene.addObject("mesh", meshdata);
    auto res = 4096;
    Containers::Array<char> data(NoInit, res * res * 4);
    Image2D tex(PixelFormat::RGBA8Unorm, {res,res}, std::move(data));
    auto view = tex.pixels<Color4ub>();
    for (int x = 0; x < res; ++x) {
        for (int y = 0; y < res; ++y) {
            auto a = x / 256;
            auto b = y / 256;
            if((a+b)%2)
                view[y][x] = Color3ub(0,0,0);
            else
                view[y][x] = Color3ub(255, 255, 255);
        }
    }

    auto plane = Primitives::grid3DSolid({1,1}, Primitives::GridFlag::GenerateTextureCoords | Primitives::GridFlag::GenerateNormals);
    auto cube = Primitives::cubeSolid();
    auto ico = Primitives::icosphereSolid(1);
    //scene.addObject("plane", upload(plane, tex));
    scene.addObject("mesh", upload(meshdata, tex));
    //scene.addNode("plane_node", "plane", ShaderType::FlatTextured);
    scene.addNode("mesh", "mesh", ShaderType::Phong);
    viewer.exec();

    return 0;

}

int main(int argc, char** argv) {
    
    ScopedTimer ti("total");

    Utility::Configuration conf(CONF_PATH, Utility::Configuration::Flag::ReadOnly);

    auto epsilon = conf.value<double>("epsilon");
    auto a = conf.value<double>("a");
    auto b = conf.value<double>("b");
    auto wMM = conf.value<double>("modica_mortola_weight");
    auto wArea = conf.value<double>("area_weight");
    auto wConnectedness = conf.value<double>("connectedness_weight");
    auto enforceConnectedness = conf.value<bool>("enforce_connectedness");
    auto postConnect = conf.value<bool>("post_connect");
    auto iterations = conf.value<int>("iterations");
    auto postConnectIterations = conf.value<int>("post_connect_iterations");
    auto inputPath = conf.value<std::string>("input_path");
    auto outputPath1 = conf.value<std::string>("output_path1");
    auto outputPath2 = conf.value<std::string>("output_path2");

    auto interfaceLoss = std::make_unique<ceres::TrivialLoss>();
    auto potentialLoss = std::make_unique<ceres::TrivialLoss>();
    auto areaLoss = std::make_unique<ceres::TrivialLoss>();
    auto connectednessLoss1 = std::make_unique<ceres::TrivialLoss>();
    auto connectednessLoss2 = std::make_unique<ceres::TrivialLoss>();

    auto meshdata = loadMeshData(inputPath);

    auto[V, F] = toEigen(meshdata);

    fmt::print("Solving Problem with {} vertices\n", V.rows());

    Eigen::VectorXd U(V.rows());
    initRandomNormal(U);

    auto* cost = new SumProblem(U.size());
    cost->push_back(std::make_unique<InterfaceEnergy>(V,F,epsilon), std::move(interfaceLoss), wMM);
    cost->push_back(std::make_unique<PotentialEnergy>(V,F,epsilon), std::move(potentialLoss), wMM);
    cost->push_back(std::make_unique<AreaRegularizer>(V,F), std::move(areaLoss), wArea);
    if(enforceConnectedness){
        cost->push_back(std::make_unique<ConnectednessConstraint>(V, F, epsilon, a, b), std::move(connectednessLoss1), wConnectedness);
        cost->push_back(std::make_unique<ConnectednessConstraint>(V, F, epsilon, -b, -a), std::move(connectednessLoss2), wConnectedness);
    }

    // Run the solver!
    ceres::GradientProblemSolver::Options options;
    options.minimizer_progress_to_stdout = true;
    options.max_num_iterations = iterations;
    options.update_state_every_iteration = true;
    ceres::GradientProblemSolver::Summary summary;

    ceres::GradientProblem problem(cost);
    ceres::Solve(options, problem, U.data(), &summary);
    std::cout << summary.BriefReport() << std::endl;
    double max = U.maxCoeff();
    double min = U.minCoeff();
    fmt::print("(min,max): ({},{})\n", min, max);
    SegmentationColorMap mapper{{{-1.05,-0.05, Color4::red()},{0.05, 1.05, Color4::blue()}}};
    writeMesh(outputPath1, V, F, U, mapper);

    if(!enforceConnectedness && postConnect){
        options.max_num_iterations = postConnectIterations;
        cost->push_back(std::make_unique<ConnectednessConstraint>(V, F, epsilon, a, b), std::move(connectednessLoss1), wConnectedness);
        cost->push_back(std::make_unique<ConnectednessConstraint>(V, F, epsilon, -b, -a), std::move(connectednessLoss2), wConnectedness);
    }

    if(postConnect){
        ceres::Solve(options, problem, U.data(), &summary);
        std::cout << summary.BriefReport() << std::endl;
    }


    writeMesh(outputPath2, V, F, U, mapper);

    auto finalCosts = cost->computeSeperateCosts(U);
    for(auto cost: finalCosts)
        std::cout << cost << std::endl;
    ScopedTimer::printStatistics();
}
