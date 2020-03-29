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

#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>

class OptimizationContext {
public:
    OptimizationContext();
    ~OptimizationContext();

    void showMenu(Magnum::ImGuiIntegration::Context&);
    void updateScene(Scene* scene);

    bool loadMesh(std::string& path);
    bool saveMesh(std::string& path);

private:

    void work();
    ceres::CallbackReturnType optimizationCallback(ceres::IterationSummary const&);

    IterationCallbackWrapper m_optimizationCallback;

    Corrade::Containers::Array<double> m_U;

    std::mutex m_meshDataMutex;
    Corrade::Containers::Optional<Trade::MeshData> m_meshData;
    bool m_reupload = false;
    bool m_newMeshData = false;

    Eigen::MatrixXi m_F;
    Eigen::MatrixXd m_V;

    double m_epsilon = 1e-2;
    double m_a = .85, m_b = .95;
    double m_weightMM = 1., m_weightConnectedness = 1., m_weightArea = 1.;

    std::string m_outputPath;
    std::string m_inputPath;

    SumProblem* m_cost = nullptr;
    ceres::GradientProblem m_problem;
    ceres::GradientProblemSolver::Options m_options;
    bool m_run = false;

    std::thread m_optThread;
    std::mutex m_mutex;
    std::condition_variable m_cv;

    std::atomic_bool m_exit = false;
};



