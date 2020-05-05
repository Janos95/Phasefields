//
// Created by janos on 04.04.20.
//

#pragma once

#include "../polygonizer/primitive_options.hpp"

#include <Magnum/Magnum.h>
#include <Magnum/Trade/MeshData.h>
#include <Corrade/Containers/Pointer.h>

#include <string>

namespace Mg = Magnum;
namespace Cr = Corrade;

enum class PrimitiveType: Magnum::UnsignedShort {
    Capsule,
    U,
    ImplicitFunction
};


bool handlePrimitive(Mg::Trade::MeshData& original, std::string& expression);


