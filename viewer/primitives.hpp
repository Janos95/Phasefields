//
// Created by janos on 04.04.20.
//

#pragma once

#include "primitive_options.hpp"

enum class PrimitiveType: Magnum::UnsignedShort {
    Capsule,
    U,
    ImplicitFunction
};

struct ComboElement {
    std::string name;
    PrimitiveType type;
    Corrade::Containers::Pointer<AbstractPrimitiveOptions> options;
};

bool handlePrimitive(Mg::Trade::MeshData& original, std::string& expression);


