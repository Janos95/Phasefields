//
// Created by janos on 20.04.20.
//

#pragma once

#include "loss_functions.hpp"
#include "shared_ressource.hpp"
#include "solver.hpp"
#include "types.hpp"
#include "meta_data.hpp"
#include "visualization_proxy.hpp"

#include <Corrade/Containers/Pointer.h>


namespace Cr = Corrade;
namespace Mg = Magnum;

template<class Scalar>
struct Functional{

    virtual bool evaluate(Scalar const* parameter, Scalar* cost, Scalar** jacobian) const = 0;
    virtual void updateInternalDataStructures() = 0;
    virtual uint32_t numParameters() const = 0;
    virtual uint32_t numResiduals() const { return 1; };
    virtual ~Functional() = default;

    Cr::Containers::Pointer<LossFunction> loss = Cr::Containers::pointer<TrivialLoss>();
    SharedRessource<Mg::Double> scaling = nullptr;
    std::string name;

    VisualizationFlags flags = {};
    VisualizationFlag* update = nullptr;

    virtual solver::Status::Value operator()(solver::IterationSummary const&) { return solver::Status::CONTINUE; };

    /*this is called from gui thread so we can update some opengl stuff if we want to */
    virtual void tickCallback() { }
    virtual void updateVisualization(VisualizationProxy& proxy) { }
    virtual void turnVisualizationOff() {}

    virtual void visualizeGradient(Cr::Containers::ArrayView<const Mg::Double> const& gradient) {}
    virtual bool drawImGuiOptions(VisualizationProxy& proxy);
};

using FunctionalF = Functional<double>;
using FunctionalD = Functional<double>;
