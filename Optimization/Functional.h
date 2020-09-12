//
// Created by janos on 20.04.20.
//

#pragma once

#include "Enums.h"
#include "Types.h"
#include "LossFunctions.h"
#include "Visualization.h"
#include "SharedRessource.h"

#include <Corrade/Containers/EnumSet.h>
#include <Corrade/Containers/Containers.h>

class adouble;

namespace Phasefield {

enum class OptionsResult : Mg::UnsignedInt {
    EvaluateProblem = 1,
    MakeUniqueVisualizer = 2,
};

using OptionsResultSet = Cr::Containers::EnumSet<OptionsResult>;

CORRADE_ENUMSET_OPERATORS(OptionsResultSet)

struct Functional {

    //enum class Traits {
    //    HasJacobian,
    //    HasHessian
    //};

    [[nodiscard]] std::size_t numParameters() const;

    void swap(Functional& other);

    explicit Functional() = default;

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

    //void
    //operator()(Containers::ArrayView<const adouble> const& params, Containers::ArrayView<const adouble> const& weights,
    //           adouble& residual) const;

    //void tapeEvaluation(double const* x) const;

    /* this is called from the gui thread so we can update some opengl stuff if we want to */

    //void updateVisualization(VisualizationProxy& proxy);

    //void turnVisualizationOff();

    OptionsResultSet drawImGuiOptions(VisualizationProxy&);

    /* mandatory */
    void* erased = nullptr;

    void (* destroy)(void*);

    size_t (* params)(void*);

    /* optional */
    //void (* vis)(void*, VisualizationProxy&) = nullptr;

    //void (* off)(void*) = nullptr;

    OptionsResultSet (* options)(void*, VisualizationProxy&) = nullptr;

    //std::size_t (* residuals)(void*) = nullptr;

    void (* evalWithGrad)(void*, ArrayView<const double> parameters,
                          ArrayView<const double> weights,
                          double& out,
                          ArrayView<double> gradP,
                          ArrayView<double> gradW);

    //void (* ad)(void*, adouble const*, adouble const*, adouble*) = nullptr;

    //int tag = -1;
    //mutable bool isFirstEvaluation = false;

    bool checkDerivatives = false;
    FunctionalType::Value functionalType = FunctionalType::Unknown;

    LossFunction loss = TrivialLoss{};
    SharedRessource<double> scaling = nullptr;

    int tag = -1;
};

}