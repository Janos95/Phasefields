//
// Created by janos on 7/1/20.
//

#pragma once

#include "Tag.h"
#include "Functional.h"
#include "Allocate.h"

#include <Corrade/Containers/ArrayView.h>
#include <new>

namespace Phasefield {

template<class T, class Lambda>
constexpr bool valid = std::is_invocable_r_v<bool, Lambda, T>;

template<class F>
Functional::Functional(F f): tag(getTag()) {

    auto checkDraw = [](auto&& f) -> decltype(f.drawImGuiOptions(std::declval<VisualizationProxy>())) {};
    auto checkSave = [](auto&& f) -> decltype(f.saveParameters(std::declval<Cr::Utility::ConfigurationGroup>())) {};
    auto checkLoad = [](auto&& f) -> decltype(f.loadParameters(std::declval<Cr::Utility::ConfigurationGroup>())) {};
    auto checkType = [](auto&& f) -> decltype(decltype(f)::type()) {};

    evalWithGrad = +[](
            void* e,
            ArrayView<const double> p,
            ArrayView<const double> w,
            double& r,
            ArrayView<double> gradP,
            ArrayView<double> gradW) {
        return (*static_cast<F*>(e))(p, w, r, gradP, gradW);
    };

#ifdef PHASEFIELD_WITH_ADOLC
    ad = +[](
            void* e,
            ArrayView<const adouble> p,
            ArrayView<const adouble> w,
            adouble& r) {
        return (*static_cast<F*>(e)).template operator()<adouble>(p, w, r, nullptr, nullptr);
    };
#endif

    /* mandatory */
    erased = allocate_buffer(sizeof(F), alignof(F));
    ::new(erased) F(std::move(f));


    destroy = +[](void* e) {
        static_cast<F*>(e)->~F();
        deallocate_buffer(e, sizeof(F), alignof(F));
    };

    params = +[](void* e) { return static_cast<F*>(e)->numParameters(); };


    /* optional */
    if constexpr(valid<F, decltype(checkDraw)>) {
        options = +[](void* e, VisualizationProxy& p) { static_cast<F*>(e)->drawImGuiOptions(p); };
    }

    if constexpr (valid<F, decltype(checkSave)>) {
        save = +[](void* e, Cr::Utility::ConfigurationGroup& group) { static_cast<F*>(e)->saveParameters(group); };
    }

    if constexpr (valid<F, decltype(checkLoad)>) {
        load = +[](void* e, Cr::Utility::ConfigurationGroup const& group) { static_cast<F*>(e)->loadParameters(group); };
    }

    if constexpr (valid<F, decltype(checkType)>)
        functionalType = F::type();

    Debug{} << tag;
}


}


