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

enum FunctionalTraits {
    HasJacobian,
    HasHessian
};

struct Functional{

    enum class Property {
        AlwaysRetape = 1,
    };

    using Properties = Cr::Containers::EnumSet<Property>;
    CORRADE_ENUMSET_FRIEND_OPERATORS(Properties)

    void updateInternalDataStructures(){
        update(erased);
        isFirstEvaluation = true;
    }

    [[nodiscard]] std::size_t numParameters() const{
        return params(erased);
    }

    [[nodiscard]] std::size_t numResiduals() const {
        if(residuals)
            return residuals(erased);
        return 1;
    };

    friend void swap(Functional& f1, Functional& f2){
        using std::swap;
        swap(f1.erased, f2.erased);
        swap(f1.destruct, f2.destruct);
        swap(f1.move, f2.move);
        swap(f1.params, f2.params);
        swap(f1.residuals, f2.residuals);
        swap(f1.vis, f2.vis);
        swap(f1.off, f2.off);
        swap(f1.update, f2.update);
        swap(f1.options, f2.options);
        swap(f1.evalWithHessian, f2.evalWithHessian);
        swap(f1.evalWithJacobian, f2.evalWithJacobian);
        swap(f1.ad, f2.ad);
        swap(f1.loss, f2.loss);
        swap(f1.scaling, f2.scaling);
        swap(f1.tag, f2.tag);
        swap(f1.isFirstEvaluation, f2.isFirstEvaluation);
        swap(f1.alwaysRetape, f2.alwaysRetape);
    }

    template<class F>
    Functional(F&& f, Property options = Property::AlwaysRetape);

    Functional(Functional const&) = delete;

    Functional& operator=(Functional const&) = delete;

    Functional(Functional&& other) noexcept{
        using std::swap;
        swap(*this, other);
    }

    Functional& operator=(Functional&& other) noexcept{
        using std::swap;
        swap(*this, other);
        return *this;
    }

    ~Functional();

    void operator()(double* params, double* residuals, SparseMatrix* jac, SparseMatrix* hess) const;

    void tapeEvaluation(double const* x) const;

    static void initTag(int& tag);

        /* this is called from the gui thread so we can update some opengl stuff if we want to */
    void updateVisualization(VisualizationProxy& proxy) {
        if(vis)
            vis(erased, proxy);
    }
    void turnVisualizationOff() {
        if(off)
            off(erased);
    }


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
