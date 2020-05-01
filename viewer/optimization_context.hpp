//
// Created by janos on 26.03.20.
//

#pragma once

#include "viewer.hpp"
#include "problem.hpp"
#include "phasefield_data.hpp"
#include "solver.hpp"
#include "connectedness_constraint.hpp"

#include <Magnum/Trade/MeshData.h>
#include <Magnum/Shaders/Flat.h>
#include <Corrade/Containers/Array.h>

#include <atomic>
#include <thread>
#include <mutex>

namespace Mg = Magnum;
namespace Cr = Corrade;

using namespace Mg::Math::Literals;

struct Path : Object3D {
    explicit Path(Object3D* parent, Trade::MeshData const& md);

    Mg::GL::Mesh mesh;
};

struct OptimizationContext final : Viewer::AbstractEventHandler {
    explicit OptimizationContext(PhasefieldData&);

    void drawImGui(Viewer&) override;
    void tickEvent(Viewer&) override;

    void drawConnectednessConstraintOptions(ConnectednessConstraint<Mg::Double>&);

    void startOptimization();
    void stopOptimization();

    solver::Status updatePhasefield(solver::IterationSummary const&);
    void appendFunctional(Cr::Containers::Array<Cr::Containers::Pointer<Functional>>&, FunctionalType);

    solver::Options options;
    std::thread thread;
    std::mutex mutex;
    std::atomic_bool optimize = false;
    std::atomic_bool update = false;

    PhasefieldData& pd;
    Cr::Containers::Array<Mg::Double> parameters;
    Cr::Containers::Array<Mg::Double> phasefield;

    SharedRessource<Mg::Double> doubleWellScaling;
    SharedRessource<Mg::Double> dirichletScaling;
    SharedRessource<Mg::Double> connectednessScaling;

    Cr::Containers::Array<Mg::Trade::MeshData> paths;
    Mg::Shaders::Flat3D pathShader;
    Mg::Color4 pathColor = Mg::Color4::green();
};


