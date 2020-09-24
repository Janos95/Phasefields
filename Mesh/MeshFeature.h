//
// Created by janos on 9/9/20.
//

#pragma once

#include "Surface.h"
#include "Types.h"

namespace Phasefield {

class MeshFeature {
public:

    explicit MeshFeature(Mesh& mesh, bool ownedByMesh = true);
    virtual ~MeshFeature();

    Mesh& mesh();

    virtual void update() = 0;

    [[nodiscard]] bool ownedByMesh() const { return m_ownedByMesh; };

    virtual char const* name() = 0;

private:
    Mesh* m_mesh = nullptr;
    bool m_ownedByMesh;
};

#define FEATURE_NAME(n) static constexpr char Name[] = n; \
char const* name() override { return Name; }

struct AngleFeature : public MeshFeature {
    using MeshFeature::MeshFeature;
    void update() override;

    FEATURE_NAME("Angle");
};

struct FaceInformationFeature : public MeshFeature {
    using MeshFeature::MeshFeature;
    void update() override;
    FEATURE_NAME("FaceInformation")
};

struct EdgeLengthFeature : public MeshFeature {
    using MeshFeature::MeshFeature;
    void update() override;
    FEATURE_NAME("EdgeLength")
};

struct GaussianCurvatureFeature : public MeshFeature {
    using MeshFeature::MeshFeature;
    void update() override;
    FEATURE_NAME("GaussianCurvature")
};

}


