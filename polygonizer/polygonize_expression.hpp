//
// Created by janos on 24.04.20.
//

#pragma once

#include <Magnum/Trade/MeshData.h>

namespace Mg = Magnum;
namespace Cr = Corrade;

struct PolygonizationOptions;

Mg::Trade::MeshData polygonizeExpression(std::string const& expression, PolygonizationOptions const& options);

