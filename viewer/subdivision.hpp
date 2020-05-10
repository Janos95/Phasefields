//
// Created by janos on 03.04.20.
//

#pragma once

#include "viewer.hpp"

void subdivide(
        std::uint32_t numSubdivisions,
        Mg::Trade::MeshData const& original,
        Cr::Containers::Array<Mg::Double>& phasefield,
        Mg::Trade::MeshData& meshData);

