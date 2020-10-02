//
// Created by janos on 8/25/20.
//

#pragma once

#include "Surface.h"
#include "Types.h"
#include "Functional.h"

class adouble;

namespace Phasefield {

struct DiffuseYamabe {

    explicit DiffuseYamabe(Mesh& m);

    template<class Scalar>
    void operator()(ArrayView<const Scalar> parameters,
                    ArrayView<const Scalar> weights,
                    Scalar& out,
                    ArrayView<Scalar> gradP,
                    ArrayView<Scalar> gradW);

    [[nodiscard]] size_t numParameters() const;

    [[nodiscard]] static FunctionalType::Value type() { return FunctionalType::DiffuseYamabe; }

    void drawImGuiOptions(VisualizationProxy& proxy);
    bool drawSolution = false;
    bool drawSolutionThresholded = false;
    bool drawSolutionGradient = false;
    double lambda = 1;

    Mesh& mesh;
};

DECLARE_FUNCTIONAL_CONSTRUCTOR(DiffuseYamabe)
DECLARE_FUNCTIONAL_OPERATOR(DiffuseYamabe, double)
DECLARE_FUNCTIONAL_OPERATOR(DiffuseYamabe, adouble)

}

