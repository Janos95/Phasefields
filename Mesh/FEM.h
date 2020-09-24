//
// Created by janos on 9/9/20.
//

#pragma once

#include "MeshFeature.h"

namespace Phasefield {

struct IntegralOperatorFeature : public MeshFeature {
    using MeshFeature::MeshFeature;

    void update() override;

    FEATURE_NAME("IntegralOperator")
};

struct MassMatrixFeature : public MeshFeature {
    using MeshFeature::MeshFeature;

    void update() override;

    FEATURE_NAME("MassMatrix")
};

struct StiffnessMatrixFeature : public MeshFeature {
    using MeshFeature::MeshFeature;

    void update() override;

    FEATURE_NAME("StiffnessMatrix")
};

struct GradientFeature : public MeshFeature {
    using MeshFeature::MeshFeature;

    void update() override;

    FEATURE_NAME("Gradient")
};

}

