//
// Created by janos on 20.04.20.
//

#pragma once

#include "types.hpp"
#include "loss_functions.hpp"

#include "optimization.h"
#include "shared_ressource.hpp"
#include "visualization_proxy.hpp"

#include <Corrade/Containers/EnumSet.h>
#include <Corrade/Containers/Containers.h>

class adouble;

enum class OptionsResult : Mg::UnsignedInt{
    EvaluateProblem = 1,
    MakeExclusive = 2,
};

using OptionsResultSet = Cr::Containers::EnumSet<OptionsResult>;
CORRADE_ENUMSET_OPERATORS(OptionsResultSet)

struct Functional{

    enum class Property {
        AlwaysRetape = 1,
    };

    using Properties = Cr::Containers::EnumSet<Property>;
    CORRADE_ENUMSET_FRIEND_OPERATORS(Properties)

    //enum class Traits {
    //    HasJacobian,
    //    HasHessian
    //};

    void updateInternalDataStructures();

    [[nodiscard]] std::size_t numParameters() const;

    friend void swap(Functional& f1, Functional& f2);

    template<class F>
    Functional(F f, Properties options = Property::AlwaysRetape);

    Functional(Functional const&) = delete;

    Functional& operator=(Functional const&) = delete;

    Functional(Functional&& other) noexcept;

    Functional& operator=(Functional&& other) noexcept;

    ~Functional();

    //void retape(Containers::StridedArrayView1D<const double> const&) const;

    void operator()(Containers::ArrayView<const double> parameters,
                    Containers::ArrayView<const double> weights,
                    double& out,
                    Containers::ArrayView<double> gradP,
                    Containers::ArrayView<double> gradW) const;

    void operator()(Containers::ArrayView<const adouble> params, Containers::ArrayView<const adouble> weights, adouble& residuals) const;

    //void tapeEvaluation(double const* x) const;

    /* this is called from the gui thread so we can update some opengl stuff if we want to */
    void updateVisualization(VisualizationProxy& proxy);

    void turnVisualizationOff();

    OptionsResultSet drawImGuiOptions();

    /* mandatory */
    void* erased = nullptr;
    void (*destroy)(void *);
    void (*move)(void*, void*);
    void (*update)(void*);
    std::size_t (*params)(void*);

    /* optional */
    void (*vis)(void*, VisualizationProxy&) = nullptr;
    void (*off)(void*) = nullptr;
    OptionsResultSet (*options)(void*) = nullptr;
    std::size_t (*residuals)(void*) = nullptr;


    void (*eval)(void*, double const*, double const*, double*, double*, double*) = nullptr;
    void (*ad)(void*, adouble const*, adouble const*, adouble*) = nullptr;
    //int tag = -1;
    mutable bool isFirstEvaluation = false;

    bool checkDerivatives = false;
    bool alwaysRetape = true;

    LossFunction loss = TrivialLoss{};
    SharedRessource<Mg::Double> scaling = nullptr;
};

