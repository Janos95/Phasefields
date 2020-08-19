//
// Created by janos on 7/1/20.
//

#pragma once

#include "Tag.h"
#include "Functional.h"

namespace Phasefield {

template<typename... Ts, typename TF>
static constexpr auto is_valid(TF) {
    return std::is_invocable<TF, Ts...>{};
}

template<class F>
Functional::Functional(F f): tag(getTag()) {

    //constexpr bool hasNumResiduals = is_valid(f)([](auto p) constexpr -> decltype(f.numResiduals()) {});

    //constexpr bool hasGradient = is_valid<double*, SparseMatrix*>(std::declval<F>())(
    //        [](auto p, auto r, auto g) constexpr -> decltype(f(p, r, g)) {});
    //constexpr bool hasHessian = is_valid<double*, double*, SparseMatrix*, SparseMatrix*>(std::declval<F>())(
    //        [](auto p, auto r, auto g, auto h) constexpr -> decltype(f(p, r, g, h)) {});

    //if constexpr(hasHessian){ /* assume we have a jacobian operator */
    //    ad = +[](void* e, adouble const* p, adouble const* w, adouble* r) {
    //        (*static_cast<F*>(e))(p, r, nullptr, nullptr);
    //    };
    //    evalWithHessian = +[](void* e, double* p, double* r, SparseMatrix* j, SparseMatrix* h) {
    //        (*static_cast<F*>(e))(p, r, j, h);
    //    };
    //} else if constexpr(hasJacobian){
    //    ad = +[](void* e, adouble* p, adouble* r) { (*static_cast<F*>(e))(p, r, nullptr); };
    //    evalWithJacobian = +[](void* e, double* p, double* r, SparseMatrix* j) { (*static_cast<F*>(e))(p, r, j); };
    //} else{
    //    ad = +[](void* e, adouble* p, adouble* r) { (*static_cast<F*>(e))(p, r); };
    //}

    //if constexpr(hasNumResiduals){
    //    residuals = +[](void* e) { return static_cast<F*>(e)->numResiduals(); };
    //}


    /* mandatory */
    erased = ::operator new(sizeof(F), std::align_val_t(alignof(F)));
    ::new(erased) F(std::move(f));


    destroy = +[](void* e) {
        static_cast<F*>(e)->~F();
        ::operator delete(e, sizeof(F), std::align_val_t(alignof(F)));
    };

    params = +[](void* e) { return static_cast<F*>(e)->numParameters(); };
    update = +[](void* e) { static_cast<F*>(e)->updateInternalDataStructures(); };


    /* optional */
    constexpr bool hasVis = is_valid(f)([](auto p) constexpr -> decltype(f.drawImGuiOptions()) {});

    if constexpr(hasVis) {
        vis = +[](void* e, VisualizationProxy& proxy) { return static_cast<F*>(e)->updateVisualization(proxy); };
        off = +[](void* e) { return static_cast<F*>(e)->turnVisualizationOff(); };
        options = +[](void* e) { return static_cast<F*>(e)->drawImGuiOptions(); };
    }

}

}

