//
// Created by janos on 20.05.20.
//

#pragma once


#include "../optimization/loss_functions.hpp"
#include "shared_ressource.hpp"
#include "../optimization/solver.hpp"
#include "types.hpp"

#include <Corrade/Containers/Pointer.h>
#include <Magnum/Trade/Trade.h>

#include <string>
#include <mutex>

struct MetaData {

    Cr::Containers::Pointer<LossFunction> loss = Cr::Containers::pointer<TrivialLoss>();
    SharedRessource<Mg::Double> scaling = nullptr;
    std::string name;

    VisualizationFlags flags = {};
    VisualizationFlag* update = nullptr;

    virtual ~MetaData() = default;
    virtual solver::Status::Value operator()(solver::IterationSummary const&) { return solver::Status::CONTINUE; };

    virtual void updateVis() { } /*this is called from gui thread so we can update some opengl stuff if we want to */
    virtual void visualizeGradient(Cr::Containers::ArrayView<const Mg::Double> const& gradient) {}
    virtual void drawImGuiOptions(bool&, DrawableType&, bool&);

    using Ptr = Cr::Containers::Pointer<MetaData>;

    template<class L>
    static Ptr fromLossAndName(L&& l, std::string n) {
        auto meta = Cr::Containers::pointer<MetaData>();
        meta->loss= Cr::Containers::pointer<std::remove_reference_t<L>>((L&&)l);
        meta->name = std::move(n);
        return meta;
    }
};


struct GradientMetaData : MetaData {

    GradientMetaData(Mg::Trade::MeshData& md, VisualizationFlags& u, std::mutex& m);

    void visualizeGradient(Cr::Containers::ArrayView<const Mg::Double> const& gradient) override;

    void drawImGuiOptions(bool& makeExclusive, DrawableType& type, bool& evaluateProblem) override;

    Mg::Trade::MeshData& meshData;
    VisualizationFlags& update;
    std::mutex& mutex;
};


struct AreaMetaData : GradientMetaData {
    Mg::Double areaRatio = 0.5;
    void drawImGuiOptions(bool&, DrawableType&, bool&) override;
};
