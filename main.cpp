

#include "ceres_cost.hpp"
#include "color_callback.hpp"
#include "phasefield_initialization.hpp"
#include "connectedness_constraint.hpp"
#include "load_mesh_data.hpp"
#include "viewer.hpp"
#include "write_mesh.hpp"

#include <scoped_timer/scoped_timer.hpp>


#include <Magnum/Trade/MeshData.h>
#include <Magnum/MeshTools/Interleave.h>
#include <Magnum/EigenIntegration/Integration.h>
#include <Magnum/Primitives/Cube.h>
#include <Magnum/Primitives/Grid.h>

#include <Corrade/Utility/Configuration.h>

#include <ceres/problem.h>
#include <ceres/solver.h>
#include <ceres/loss_function.h>
#include <ceres/gradient_problem.h>
#include <ceres/gradient_problem_solver.h>

#include <thread>
#include <algorithm>
#include <iostream>

using namespace Magnum;
using namespace Corrade;
using namespace Corrade::Containers;

using namespace std::chrono_literals;



int main(int argc, char** argv) {
    
    std::cout << CONF_PATH << std::endl;

    ScopedTimer ti("total");
    
    Utility::Configuration conf(CONF_PATH, Utility::Configuration::Flag::ReadOnly);
    auto eps = conf.value<double>("epsilon");
    std::cout << eps << std::endl;
    return 1;

    constexpr double epsilon = 1e-2;

    auto interfaceLoss = std::make_unique<ceres::TrivialLoss>();
    auto potentialLoss = std::make_unique<ceres::TrivialLoss>();
    auto areaLoss = std::make_unique<ceres::TrivialLoss>();
    auto connectednessLoss1 = std::make_unique<ceres::TrivialLoss>();
    auto connectednessLoss2 = std::make_unique<ceres::TrivialLoss>();

    auto meshdata = loadMeshData("/home/janos/data/centaur.ply");
    //auto meshdata = Primitives::grid3DSolid({20,20});

    auto varr = meshdata.positions3DAsArray();
    auto farr = meshdata.indicesAsArray();

    printf("Solving Problem with %d vertices\n", varr.size());
    //Debug{} << varr;
    //Debug{} << farr;

    //Debug{} << varr.size();
    //Debug{} << farr.size();

    using MatrixXU = Eigen::Matrix<UnsignedInt, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    CORRADE_INTERNAL_ASSERT(farr.size() % 3 == 0);
    Eigen::MatrixXi F = Eigen::Map<const MatrixXU>(farr.data(), farr.size() / 3, 3).cast<int>();

    Eigen::MatrixXd V(varr.size(), 3);
    for (int i = 0; i < meshdata.vertexCount(); ++i) {
       V.row(i) = EigenIntegration::cast<Eigen::Vector3f>(varr[i]).cast<double>();
    }

    //normalize vertex length to be in [0,1] for better numerical stability
    double maxVertexNorm = V.rowwise().norm().maxCoeff();
    V *= 1. / maxVertexNorm;

    //std::cout << V << std::endl;
    //std::cout << F << std::endl;

    //std::cout << V.rows() << std::endl;
    //std::cout << F.rows() << std::endl;

    Eigen::VectorXd U(V.rows());
    initRandomNormal(U);

    auto* cost = new SumProblem(U.size());
    cost->push_back(std::make_unique<InterfaceEnergy>(V,F,epsilon), std::move(interfaceLoss));
    cost->push_back(std::make_unique<PotentialEnergy>(V,F,epsilon), std::move(potentialLoss));
    cost->push_back(std::make_unique<AreaRegularizer>(V,F), std::move(areaLoss));
    //cost->push_back(std::make_unique<ConnectednessConstraint>(V, F, epsilon, 0.85, 0.95), std::move(connectednessLoss1));
    //cost->push_back(std::make_unique<ConnectednessConstraint>(V, F, epsilon, -0.95, -0.85), std::move(connectednessLoss2));

    //Debug{} << varr;
    //Debug{} << farr;

    //Viewer viewer(argc,argv);
    //Debug{} << meshdata.vertexData().size();
    //viewer.scene.addObject("mesh", meshdata, CompileFlag::AddColorAttribute|CompileFlag::GenerateFlatNormals);
    //auto& vertices = viewer.scene.getObject("mesh")->vertices;
    //Containers::Array<char> data(vertices.size());
    //Utility::copy({vertices.mapRead(), vertices.size()}, data);
    //ColorCallback colorCB(U, std::move(data));

    //{
    //    colorCB(viewer.scene);
    //    data = Containers::Array<char>(vertices.size());
    //    Containers::ArrayView raw(vertices.mapRead(), vertices.size());
    //    CORRADE_INTERNAL_ASSERT(raw);
    //    Containers::StridedArrayView1D<void> erasedViewC(
    //            data,
    //            data.begin() + 6 * sizeof(Float) /*start*/,
    //            U.size() /*size */ ,
    //            10 * sizeof(Float) /* stride */);
    //    auto viewC = arrayCast<Color4>(erasedViewC);

    //    Containers::StridedArrayView1D<void> erasedViewP(
    //            data,
    //            data.begin()  /*start*/,
    //            U.size() /*size */ ,
    //            10 * sizeof(Float) /* stride */);
    //    auto viewP = arrayCast<Vector3>(erasedViewP);

    //    Containers::StridedArrayView1D<void> erasedViewN(
    //            data,
    //            data.begin() + 3 * sizeof(Float) /*start*/,
    //            U.size() /*size */ ,
    //            10 * sizeof(Float) /* stride */);
    //    auto viewN = arrayCast<Vector3>(erasedViewN);

    //    Debug{} << viewP;
    //    Debug{} << viewN;
    //    Debug{} << viewC;
    //    vertices.unmap();
    //}

    // Run the solver!
    ceres::GradientProblemSolver::Options options;
    options.minimizer_progress_to_stdout = true;
    options.max_num_iterations = 10000;
    options.update_state_every_iteration = true;
    ceres::GradientProblemSolver::Summary summary;

    ceres::GradientProblem problem(cost);
    ceres::Solve(options, problem, U.data(), &summary);
    std::cout << summary.BriefReport() << std::endl;

    writeMesh("/home/janos/data/centaur_out1.ply", varr, farr, U, true);
    //std::thread t([&]{
    //    ceres::Solve(options, problem, U.data(), &summary);
    //    std::cout << summary.BriefReport() << '\n';
    //});

    //std::thread t([&]{
    //    while(1){
    //        initRandomNormal(U);
    //        colorCB(ceres::IterationSummary{});
    //        std::this_thread::sleep_for(1s/60);
    //    }
    //});

    // Plot the mesh
    //viewer.callbacks.emplace_back(colorCB);
    //viewer.exec();

    //t.join();


    auto finalCosts = cost->computeSeperateCosts(U);
    for(auto cost: finalCosts)
        std::cout << cost << std::endl;
    ScopedTimer::printStatistics();
}
