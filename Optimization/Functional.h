//
// Created by janos on 20.04.20.
//

#pragma once

#include "Enums.h"
#include "Types.h"
#include "LossFunctions.h"
#include "Visualization.h"

#include <Corrade/Utility/Utility.h>

class adouble;

namespace Phasefield {

struct Functional {

    //enum class Traits {
    //    HasJacobian,
    //    HasHessian
    //};

    [[nodiscard]] std::size_t numParameters() const;

    void swap(Functional& other);

    explicit Functional() = default;
    explicit Functional(std::nullptr_t) {}

    template<class F>
    /* implicit */ Functional(F f);

    Functional(Functional const&) = delete;

    Functional& operator=(Functional const&) = delete;

    Functional(Functional&& other) noexcept;

    Functional& operator=(Functional&& other) noexcept;

    ~Functional();

    //void retape(Containers::StridedArrayView1D<const double> const&) const;

    void operator()(ArrayView<const double> parameters,
                    ArrayView<const double> weights,
                    double& out,
                    ArrayView<double> gradP,
                    ArrayView<double> gradW) const;

    void loadParameters(Cr::Utility::ConfigurationGroup const&);

    void saveParameters(Cr::Utility::ConfigurationGroup&) const;

    void drawImGuiOptions(VisualizationProxy&);

    /* mandatory */
    void* erased = nullptr;

    void (* destroy)(void*);

    size_t (* params)(void*);

    /* optional */
    void (* load)(void*, Cr::Utility::ConfigurationGroup const&) = nullptr;
    void (* save)(void*, Cr::Utility::ConfigurationGroup&) = nullptr;

    //void (* off)(void*) = nullptr;

    void (* options)(void*, VisualizationProxy&) = nullptr;

    //std::size_t (* residuals)(void*) = nullptr;

    void (* evalWithGrad)(void*, ArrayView<const double> parameters,
                          ArrayView<const double> weights,
                          double& out,
                          ArrayView<double> gradP,
                          ArrayView<double> gradW);

    void (* ad)(void*, ArrayView<const adouble>, ArrayView<const adouble>, adouble&) = nullptr;

    bool check(ArrayView<const double>) const;

    //int tag = -1;
    //mutable bool isFirstEvaluation = false;

    bool checkDerivatives = false;
    FunctionalType::Value functionalType = FunctionalType::Unknown;

    LossFunction loss = TrivialLoss{};
    double* scaling = nullptr;

    size_t tag = Invalid;
    bool drawGradient = false;
    bool disable = false;
};


#define DECLARE_FUNCTIONAL_OPERATOR(name, type) \
extern template void name::operator()(ArrayView<const type> parameters, \
                                      ArrayView<const type> weights, \
                                      type& out, \
                                      ArrayView<type> gradP, \
                                      ArrayView<type> gradW);

#define DEFINE_FUNCTIONAL_OPERATOR(name, type) \
template void name::operator()(ArrayView<const type> parameters, \
                                      ArrayView<const type> weights, \
                                      type& out, \
                                      ArrayView<type> gradP, \
                                      ArrayView<type> gradW);

#define DECLARE_FUNCTIONAL_CONSTRUCTOR(name) extern template Functional::Functional(name);
#define DEFINE_FUNCTIONAL_CONSTRUCTOR(name) template Functional::Functional(name);

}

