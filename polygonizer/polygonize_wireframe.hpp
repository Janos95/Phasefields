//
// Created by janos on 08.05.20.
//

#pragma once

#include "primitive_options.hpp"

#include <Corrade/Containers/ArrayView.h>
#include <Magnum/Trade/MeshData.h>

namespace Mg = Magnum;
namespace Cr = Corrade;

Mg::Trade::MeshData polygonizeWireframe(
        Cr::Containers::ArrayView<const Mg::Vector3d> const& vertices,
        Cr::Containers::ArrayView<const Mg::Vector3ui> const& triangles,
        PolygonizationOptions const& options);


