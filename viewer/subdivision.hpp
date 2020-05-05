//
// Created by janos on 03.04.20.
//

#pragma once

#include "viewer.hpp"

void subdivide(
        std::uint32_t numSubdivisions,
        Trade::MeshData& data,
        Containers::Array<Double>& phasefield,
        Containers::Array<Vector3d>& V,
        Containers::Array<UnsignedInt>& T,
        Trade::MeshData& meshData);

