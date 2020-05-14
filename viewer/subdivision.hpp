//
// Created by janos on 03.04.20.
//

#pragma once

#include <Magnum/Magnum.h>
#include <Corrade/Containers/Containers.h>

namespace Cr = Corrade;
namespace Mg = Magnum;

void subdivide(
        std::uint32_t numSubdivisions,
        Cr::Containers::Array<Mg::UnsignedInt>& indices,
        Cr::Containers::Array<Mg::Vector3d>& vertices,
        Cr::Containers::Array<Mg::Double>& phasefield);