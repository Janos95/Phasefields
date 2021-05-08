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

template<class F>
Functional::Functional(F f): tag(getTag()) {

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
    if constexpr(requires(VisualizationProxy& p) { f.drawImGuiOptions(p); }) {
        options = +[](void* e, VisualizationProxy& p) { static_cast<F*>(e)->drawImGuiOptions(p); };
    }

    if constexpr (requires(Cr::Utility::ConfigurationGroup& group) { f.saveParameters(group); }) {
        save = +[](void* e, Cr::Utility::ConfigurationGroup& group) { static_cast<F*>(e)->saveParameters(group); };
    }

    if constexpr (requires(Cr::Utility::ConfigurationGroup const& group) { f.loadParameters(group); }) {
        load = +[](void* e, Cr::Utility::ConfigurationGroup const& group) { static_cast<F*>(e)->loadParameters(group); };
    }

    if constexpr (requires { F::type(); })
        functionalType = F::type();

    Debug{} << tag;
}


}


