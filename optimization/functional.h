//
// Created by janos on 20.04.20.
//

#pragma once

#include "types.hpp"
#include "loss_functions.hpp"

#include "optimization.h"
#include "shared_ressource.hpp"
#include "../visualization_proxy/visualization_proxy.hpp"

#include <Corrade/Containers/EnumSet.h>

class adouble;

enum class ImGuiResult : Mg::UnsignedInt{
    EvaluateProblem = 1,
    MakeExclusive = 2,
};
using ImGuiResults = Cr::Containers::EnumSet<ImGuiResult>;
CORRADE_ENUMSET_OPERATORS(ImGuiResults)


struct Functional{

    enum class Property {
        AlwaysRetape = 1,
    };

    using Properties = Cr::Containers::EnumSet<Property>;
    CORRADE_ENUMSET_FRIEND_OPERATORS(Properties)

    enum class Traits {
        HasJacobian,
        HasHessian
    };


    void updateInternalDataStructures();

    [[nodiscard]] std::size_t numParameters() const;

    [[nodiscard]] std::size_t numResiduals() const;

    friend void swap(Functional& f1, Functional& f2);

    template<class F>
    Functional(F f, Property options = Property::AlwaysRetape);

    Functional(Functional const&) = delete;

    Functional& operator=(Functional const&) = delete;

    Functional(Functional&& other) noexcept;

    Functional& operator=(Functional&& other) noexcept;

    ~Functional();

    void operator()(double* params, double* residuals, SparseMatrix* jac, SparseMatrix* hess) const;

    void tapeEvaluation(double const* x) const;

    static void initTag(int& tag);

    /* this is called from the gui thread so we can update some opengl stuff if we want to */
    void updateVisualization(VisualizationProxy& proxy);

    void turnVisualizationOff();

    ImGuiResults drawImGuiOptions();

    /* mandatory */
    void* erased = nullptr;
    void (*destruct)(void *);
    void (*move)(void*, void*);
    void (*update)(void*);
    std::size_t (*params)(void*);

    /* optional */
    void (*vis)(void*, VisualizationProxy&) = nullptr;
    void (*off)(void*) = nullptr;
    bool (*options)(void*) = nullptr;
    std::size_t (*residuals)(void*) = nullptr;


    void (*evalWithJacobian)(void*, double*, double*, SparseMatrix*) = nullptr;
    void (*evalWithHessian)(void*, double*, double*, SparseMatrix*, SparseMatrix*) = nullptr;
    void (*ad)(void*, adouble*, adouble*) = nullptr;
    int tag = -1;
    mutable bool isFirstEvaluation = false;

    bool useAd = false;
    bool alwaysRetape = true;

    LossFunction loss = TrivialLoss{};
    SharedRessource<Mg::Double> scaling = nullptr;
};

