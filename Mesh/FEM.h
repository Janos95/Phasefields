//
// Created by janos on 9/9/20.
//

#pragma once

#include "MeshFeature.h"

namespace Phasefield {

struct IntegralOperatorFeature : public MeshFeature {
    using MeshFeature::MeshFeature;

    void update() override;

    const char* name() override { return "IntegralOperator"; }
};

struct MassMatrixFeature : public MeshFeature {
    using MeshFeature::MeshFeature;

    void update() override;

    const char* name() override { return "MassMatrix"; }
};

//struct StiffnessMatrixFeature : public MeshFeature {
//    using MeshFeature::MeshFeature;
//
//    void update() override;
//
//    const char* name() override { return "StiffnessMatrix"; }
//};

struct GradientFeature : public MeshFeature {
    using MeshFeature::MeshFeature;

    void update() override;

    const char* name() override { return "Gradient"; }
};

}

