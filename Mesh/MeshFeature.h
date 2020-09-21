//
// Created by janos on 9/9/20.
//

#pragma once

#include "Surface.h"
#include "Types.h"

namespace Phasefield {

class MeshFeature {
public:

    MeshFeature(Mesh& mesh);

    Mesh& mesh();

    virtual void update() = 0;

    virtual const char* name() = 0;

private:
    Mesh& m_mesh;
};

struct AngleFeature : public MeshFeature {
    using MeshFeature::MeshFeature;
    void update() override;
    const char* name() override { return "Angle"; }
};

struct FaceInformationFeature : public MeshFeature {
    using MeshFeature::MeshFeature;
    void update() override;
    const char* name() override { return "FaceArea"; }
};

struct EdgeLengthFeature : public MeshFeature {
    using MeshFeature::MeshFeature;
    void update() override;
    const char* name() override { return "EdgeLength"; }
};

struct GaussianCurvatureFeature : public MeshFeature {
    using MeshFeature::MeshFeature;
    void update() override;
    const char* name() override { return "GaussianCurvature"; }
};

}


