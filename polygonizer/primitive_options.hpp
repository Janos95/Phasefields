//
// Created by janos on 29.04.20.
//

#pragma once

#include <cstdint>

struct AbstractPrimitiveOptions { virtual ~AbstractPrimitiveOptions() = default; };

struct CapsuleOptions : AbstractPrimitiveOptions {
    std::uint32_t hemisphereRings = 10; //	Number of (face) rings for each hemisphere. Must be larger or equal to 1.
    std::uint32_t cylinderRings = 30; // 	Number of (face) rings for cylinder. Must be larger or equal to 1.
    std::uint32_t segments = 30; //	Number of (face) segments. Must be larger or equal to 3.
    float radius = 1.f;
    float length = 5.f;
};

struct UOptions : AbstractPrimitiveOptions {
    float height = 1.f;
    float width = 1.f;
    float innerWidth = .5f;
    float innerHeight = .5f;
};

struct PolygonizationOptions : AbstractPrimitiveOptions {
    float angleBound = 30.f;
    float radiusBound = 0.1f;
    float distanceBound = 0.1f;
    float boundingSphereRadius = 2.f;
};
