//
// Created by janos on 26.03.20.
//

#pragma once

#include "viewer.hpp"
#include "problem.hpp"
#include "ray_mesh_intersection.hpp"

#include <Magnum/ImGuiIntegration/Context.h>

#include <ceres/gradient_problem.h>
#include <ceres/gradient_problem_solver.h>
#include <ceres/iteration_callback.h>

#include <Eigen/Core>

#include <mutex>
#include <thread>

class OptimizationContext : public Viewer::AbstractEventHandler {
public:
    OptimizationContext();
    ~OptimizationContext();

    void drawEvent() override;

    void startOptimization();
    void stopOptimization();

    enum class InitializationFlag: Magnum::UnsignedByte {
        RandomNormal = 1
    };

    void initializePhasefield(InitializationFlag flag);

private:

    ceres::CallbackReturnType optimizationCallback(ceres::IterationSummary const&);

    IterationCallbackWrapper m_optimizationCallback;

    Eigen::MatrixXi m_F;
    Eigen::MatrixXd m_V;

    double m_epsilon = .05;
    float m_posA = .85f, m_posB = .95f;
    float m_negA = -.95f, m_negB = -.85f;

    SumProblem* m_cost = nullptr;
    Corrade::Containers::Optional<ceres::GradientProblem> m_problem;
    ceres::GradientProblemSolver::Options m_options;

    std::vector<std::pair<std::string, double>> m_weights =
            {
                {"Modica Mortola", 1.},
                {"Area Regularization", 1.},
                {"Connectedness Constraint Positive", 1.},
                {"Connectedness Constraint Negative", 1.}
            };
};



