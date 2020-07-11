//
// Created by janos on 7/1/20.
//

#pragma once

#include "functional.h"


template <typename... Ts, typename TF>
static constexpr auto is_valid(TF){
    return std::is_invocable<TF, Ts...>{};
}

template<class F>
Functional::Functional(F&& f, Options options){
    using T = std::remove_reference_t<F>;

    erased = std::malloc(sizeof(T));
    ::new(erased) F((F&&)f);

    constexpr bool hasNumResiduals = is_valid(f)([](auto p) constexpr -> decltype(f.numResiduals()) {});
    constexpr bool hasVis = is_valid(f)([](auto p) constexpr -> decltype(f.drawImGuiOptions()) {});

    constexpr bool hasJacobian = is_valid<double*, SparseMatrix*>(std::declval<T>())([](auto p, auto r, auto g) constexpr -> decltype(f(p,r,g)) {});
    constexpr bool hasHessian = is_valid<double*, double*, SparseMatrix*, SparseMatrix*>(std::declval<T>())([](auto p, auto r, auto g, auto h) constexpr -> decltype(f(p,r,g,h)) {});

    if constexpr(hasHessian) { /* assume we have a jacobian operator */
        ad = +[](void *e, adouble *p, adouble *r) { (*static_cast<T*>(e))(p, r, nullptr, nullptr); };
        evalWithHessian = +[](void *e, double *p, double *r, SparseMatrix* j, SparseMatrix* h) { (*static_cast<T *>(e))(p, r, j, h); };
    } else if constexpr(hasJacobian) {
        ad = +[](void *e, adouble *p, adouble *r) { (*static_cast<T*>(e))(p, r, nullptr); };
        evalWithJacobian = +[](void *e, double *p, double *r, SparseMatrix* j) { (*static_cast<T *>(e))(p, r, j); };
    } else {
        ad = +[](void *e, adouble *p, adouble *r) { (*static_cast<T*>(e))(p, r); };
    }

    initTag(tag);

    if constexpr(hasNumResiduals){
        residuals = +[](void* e) { return static_cast<T*>(e)->numResiduals(); };
    }

    if constexpr(hasVis){
        residuals = +[](void* e) { return static_cast<T*>(e)->numResiduals(); };
        vis = +[](void* e, VisualizationProxy& proxy) { return static_cast<T*>(e)->updateVisualization(proxy); };
        off = +[](void* e) { return static_cast<T*>(e)->turnVisualizationOff(); };
        options = +[](void* e) { return static_cast<T*>(e)->drawImGuiOptions(); };
    }

    destruct = +[](void* e) { return static_cast<T*>(e)->~T(); };
}