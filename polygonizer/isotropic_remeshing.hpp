//
// Created by janos on 12.05.20.
//

#pragma once

#include <Corrade/Containers/Containers.h>
#include <Magnum/Magnum.h>

namespace Mg = Magnum;
namespace Cr = Corrade;

void isotropicRemeshing(
        Cr::Containers::Array<Mg::Vector3d>& vertices,
        Cr::Containers::Array<Mg::UnsignedInt>& indices);


