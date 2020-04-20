//
// Created by janos on 26.03.20.
//

#pragma once

#include "viewer.hpp"
#include "problem.hpp"
#include "phasefield_data.hpp"
#include "solver.hpp"

#include <Corrade/Containers/Array.h>

#include <tbb/task_group.h>

#include <atomic>

namespace Mg = Magnum;
namespace Cr = Corrade;

struct OptimizationContext final : Viewer::AbstractEventHandler {
    explicit OptimizationContext(PhasefieldData&);

    void drawImGui(Viewer&) override;
    void tickEvent(Viewer&) override;

    void startOptimization();
    void stopOptimization();

    Problem problem;
    solver::Options options;
    Corrade::Containers::Array<Mg::Double> phasefield;
    tbb::task_group taskGroup;

    std::atomic_bool optimize = false;

    PhasefieldData& pd;

    bool showShortestPaths = false;
};



