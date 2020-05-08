//
// Created by janos on 03.04.20.
//

#pragma once

#include "viewer.hpp"

void subdivide(
        std::uint32_t numSubdivisions,
        Trade::MeshData const& original,
        Containers::Array<Double>& phasefield,
        Trade::MeshData& meshData);

