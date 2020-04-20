//
// Created by janos on 28.03.20.
//

#pragma once

#include "unique_function.h"
#include "functional.hpp"

#include <Corrade/Containers/GrowableArray.h>

namespace Cr = Corrade;

//struct IterationCallbackWrapper : public ceres::IterationCallback{
//
//    template<class F>
//    IterationCallbackWrapper(F&& f) : callback((F&&)f){}
//
//    ceres::CallbackReturnType operator()(ceres::IterationSummary const& summary) override {
//         return callback(summary);
//    }
//
//     unique_function<ceres::CallbackReturnType(ceres::IterationSummary const&)> callback;
//};
namespace solver {

struct Problem {

    struct Cost {
        std::string_view name;
        double cost;
    };

    [[nodiscard]] Cr::Containers::Array<Cost> evaluateFunctionals(double const *parameters) const;

    bool evaluate(const double *parameters,
                  double *cost,
                  double *jacobians) const;

    [[nodiscard]] int numParameters() const;

    Cr::Containers::Array<Cr::Containers::Pointer<Functional>> functionals;
    Cr::Containers::Array<Cr::Containers::Pointer<Functional>> constraints;
};

}