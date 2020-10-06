//
// Created by janos on 8/25/20.
//

#pragma once

#include "Surface.h"
#include "Types.h"
#include "Functional.h"

class adouble;

namespace Phasefield {

SMART_ENUM(EnergyType, size_t, Dirichlet, Hencky);

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
    bool curvatureRescaling = false;
    double* lambdaScaling = nullptr;
    double lambdaWeight = 1;

    Mesh& mesh;


    double getRescalingFactor(Face f) const;
    double getRescalingFactor(Vertex v) const;

    EnergyType::Value energy = EnergyType::Hencky;
};

DECLARE_FUNCTIONAL_CONSTRUCTOR(DiffuseYamabe)
DECLARE_FUNCTIONAL_OPERATOR(DiffuseYamabe, double)
DECLARE_FUNCTIONAL_OPERATOR(DiffuseYamabe, adouble)

}
