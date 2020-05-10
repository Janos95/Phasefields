//
// Created by janos on 20.04.20.
//

#pragma once

#include "loss_functions.hpp"
#include "shared_ressource.hpp"
#include "solver.hpp"
#include "types.hpp"

#include <Corrade/Containers/Containers.h>
#include <mutex>

namespace Cr = Corrade;
namespace Mg = Magnum;

enum class FunctionalType : Mg::UnsignedInt {
    Undefined = 0,
    DirichletEnergy = 1,
    DoubleWellPotential = 2,
    Area1 = 3,
    Connectedness = 4,
    Area2 = 5,
};

struct Functional{

    struct MetaData {

        Cr::Containers::Pointer<LossFunction> loss = Cr::Containers::pointer<TrivialLoss>();
        SharedRessource<Mg::Double> scaling = nullptr;

        VisualizationFlags flags = {};

        VisualizationFlag* update;
        std::mutex* mutex;

        virtual ~MetaData() = default;
        virtual solver::Status operator()(solver::IterationSummary const&) { return solver::Status::CONTINUE; };

        virtual void updateVis() { } /*this is called from gui thread so we can update some opengl stuff if we want to */
        virtual void visualizeGradient(Cr::Containers::ArrayView<const Mg::Double> const& gradient) {};

        using Ptr = Cr::Containers::Pointer<MetaData>;

        template<class L>
        static Ptr AllocateFromLoss(L&& l) {
            auto meta = Cr::Containers::pointer<MetaData>();
            meta->loss= Cr::Containers::pointer<std::remove_reference_t<L>>((L&&)l);
            return meta;
        }
    };

    Functional(Cr::Containers::Pointer<MetaData> meta, FunctionalType t) : metaData(std::move(meta)), type(t)
    {
    }

    mutable Cr::Containers::Pointer<MetaData> metaData = Cr::Containers::pointer<MetaData>();
    FunctionalType type;

    virtual bool evaluate(double const* parameter, double* cost, double* jacobian) const = 0;
    virtual int numParameters() const = 0;
    virtual ~Functional() = default;
};



