//
// Created by janos on 8/25/20.
//

#pragma once

#include "Mesh.h"
#include "Surface.h"
#include "Types.h"
#include "Functional.h"

#include <Corrade/Containers/Array.h>

class adouble;

namespace Phasefield {

struct LinearizedElasticity {

    explicit LinearizedElasticity(Mesh& mesh);

    template<class Scalar>
    void operator()(ArrayView<const Scalar> parameters,
                    [[maybe_unused]] ArrayView<const Scalar> weights,
                    Scalar& out,
                    ArrayView<Scalar> gradP,
                    [[maybe_unused]] ArrayView<Scalar> gradW);

    [[nodiscard]] size_t numParameters() const;

    [[nodiscard]] static FunctionalType::Value type() { return FunctionalType::DiffuseYamabe; }

    void drawImGuiOptions(VisualizationProxy& proxy);

    double surfaceLoad();
    void computeDisplacements();

    double youngModulus;
    double poissonRatio;
    double forceNorm;

    VertexData<double> displacementsX, displacementsY;
    VertexData<Vector2d> positions;

    Array<Vertex> neumannBoundary;
    Array<Vertex> dirichletBoundary;

    Array<double> neumannElements;
    Array<double> outwardNormals;

    bool applyDisplacements;

    Mesh& mesh;
};

DECLARE_FUNCTIONAL_CONSTRUCTOR(LinearizedElasticity)
DECLARE_FUNCTIONAL_OPERATOR(LinearizedElasticity, double)

#ifdef PHASEFIELD_WITH_ADOLC
DECLARE_FUNCTIONAL_OPERATOR(LinearizedElasticity, adouble)
#endif

}

