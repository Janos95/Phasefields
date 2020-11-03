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

    explicit DiffuseYamabe(Mesh& mesh);

    template<class Scalar>
    void operator()(ArrayView<const Scalar> parameters,
                    [[maybe_unused]] ArrayView<const Scalar> weights,
                    Scalar& out,
                    ArrayView<Scalar> gradP,
                    [[maybe_unused]] ArrayView<Scalar> gradW);

    [[nodiscard]] size_t numParameters() const;

    [[nodiscard]] static FunctionalType::Value type() { return FunctionalType::DiffuseYamabe; }

    void drawImGuiOptions(VisualizationProxy& proxy);
    bool drawSolution = false;
    bool drawSolutionThresholded = false;
    bool drawSolutionGradient = false;
    bool curvatureRescaling = false;
    double lambdaWeight = 1;
    double* scaling = nullptr;

    Mesh& mesh;

    double getRescalingFactor(Face f) const;
    double getRescalingFactor(Vertex v) const;

    EnergyType::Value energy = EnergyType::Hencky;
    bool positivePhase = true;
    bool negativePhase = true;
};

DECLARE_FUNCTIONAL_CONSTRUCTOR(DiffuseYamabe)
DECLARE_FUNCTIONAL_OPERATOR(DiffuseYamabe, double)

#ifdef PHASEFIELD_WITH_ADOLC
DECLARE_FUNCTIONAL_OPERATOR(DiffuseYamabe, adouble)
#endif

}

