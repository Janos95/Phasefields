//
// Created by janos on 26.03.20.
//

#pragma once

#include "viewer.hpp"
#include "problem.hpp"

#include <Magnum/ImGuiIntegration/Context.h>

#include <ceres/gradient_problem.h>
#include <ceres/gradient_problem_solver.h>
#include <ceres/iteration_callback.h>

#include <Eigen/Core>

#include <folly/futures/Future.h>
#include <folly/executors/ThreadedExecutor.h>

#include <mutex>
#include <thread>


class OptimizationContext {
public:
    OptimizationContext();
    ~OptimizationContext();

    void showMenu(Magnum::ImGuiIntegration::Context&);
    void updateScene(Scene*& scene);

    bool loadMesh(std::string const& path);
    bool saveMesh(std::string const& path);

    void startOptimization();
    void stopOptimization();

    enum class InitializationFlag: Magnum::UnsignedByte {
        RandomNormal = 1
    };

    void initializePhasefield(InitializationFlag flag);

private:

    ceres::CallbackReturnType optimizationCallback(ceres::IterationSummary const&);

    IterationCallbackWrapper m_optimizationCallback;

    Eigen::VectorXd m_U;

    folly::Future<folly::Unit> m_optimizationFuture;
    folly::ThreadedExecutor m_executor;

    std::mutex m_mutex;
    Corrade::Containers::Optional<Trade::MeshData> m_meshData;
    std::atomic_bool m_reupload = false;
    std::atomic_bool m_newMeshData = false;
    std::atomic_bool m_exit = false;

    Eigen::MatrixXi m_F;
    Eigen::MatrixXd m_V;

    double m_epsilon = 1e-2;
    float m_posA = .85f, m_posB = .95f;
    float m_negA = -.95f, m_negB = -.85f;

    std::vector<char> m_outputPath;
    std::vector<char> m_inputPath;

    SumProblem* m_cost = nullptr;
    Corrade::Containers::Optional<ceres::GradientProblem> m_problem;
    ceres::GradientProblemSolver::Options m_options;

    std::vector<std::pair<std::string, double>> m_checkBoxes = {{"Modica Mortola", 1.},{"Area Regularization", 1.},{"Connectedness Constraint Positive", 1.},{"Connectedness Constraint Negative", 1.}};
};



