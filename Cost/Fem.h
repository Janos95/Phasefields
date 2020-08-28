//
// Created by janos on 30.04.20.
//

#pragma once

#include "SparseMatrix.h"

#include <Corrade/Containers/Containers.h>
#include <Corrade/Containers/Array.h>
#include <Magnum/Magnum.h>
#include <Magnum/Math/Math.h>

namespace Phasefield {

namespace Cr = Corrade;
namespace Mg = Magnum;

Mg::Double computeArea(Mg::Vector3d const& a, Mg::Vector3d const& b, Mg::Vector3d const& c);

Cr::Containers::Array<Mg::Double> computeAreas(
        const Cr::Containers::ArrayView<const Mg::Vector3ui>& triangles,
        const Cr::Containers::ArrayView<const Mg::Vector3d>& vertices);

Cr::Containers::Array<Mg::Double> computeIntegralOperator(
        const Cr::Containers::ArrayView<const Mg::Vector3ui>& triangles,
        const Cr::Containers::ArrayView<const Mg::Vector3d>& vertices);

Cr::Containers::Array<Triplet> computeMassMatrix(
        const Cr::Containers::ArrayView<const Mg::Vector3ui>& triangles,
        const Cr::Containers::ArrayView<const Mg::Vector3d>& vertices);

Cr::Containers::Array<Mg::Vector3d> gradient(
        const Cr::Containers::ArrayView<const Mg::Vector3ui>& triangles,
        const Cr::Containers::ArrayView<const Mg::Vector3d>& vertices);

Cr::Containers::Array<Triplet> computeStiffnessMatrix(
        const Cr::Containers::ArrayView<const Mg::Vector3ui>& triangles,
        const Cr::Containers::ArrayView<const Mg::Vector3d>& vertices);
}