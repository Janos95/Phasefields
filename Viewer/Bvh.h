//
// Created by janos on 10/11/20.
//

#pragma once

#include "Surface.h"
#include "MeshFeature.h"

namespace Phasefield {

struct BVHAdapter : MeshFeature {

    struct Impl;
    Impl* impl;

    explicit BVHAdapter(Mesh& mesh);
    ~BVHAdapter();

    struct Intersection {
        float t, u, v;
        size_t idx;
    };

    bool computeIntersection(Vector3 const& p, Vector3 const& dir, Intersection& intersection);

    void update() override;

    FEATURE_NAME("BVH")

};

}

