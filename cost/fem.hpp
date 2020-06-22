//
// Created by janos on 30.04.20.
//

#pragma once

#include <Eigen/SparseCore>

#include <Corrade/Containers/Containers.h>
#include <Magnum/Magnum.h>
#include <Magnum/Math/Math.h>

namespace Cr = Corrade;
namespace Mg = Magnum;

Mg::Double computeArea(Mg::Vector3d const& a, Mg::Vector3d const& b, Mg::Vector3d const& c);

Cr::Containers::Array<Mg::Double> computeAreas(
        const Cr::Containers::ArrayView<const Mg::Vector3ui>& indices,
        const Cr::Containers::ArrayView<const Mg::Vector3d>& data);

Cr::Containers::Array<Mg::Double> computeIntegralOperator(
        const Cr::Containers::ArrayView<const Mg::Vector3ui>& indices,
        const Cr::Containers::ArrayView<const Mg::Vector3d>& data);

Eigen::SparseMatrix<Mg::Double> computeMassMatrix(
        const Cr::Containers::ArrayView<const Mg::Vector3ui>& triangles,
        const Cr::Containers::ArrayView<const Mg::Vector3d>& vertices);

Eigen::SparseMatrix<Mg::Double> gradient(
        const Cr::Containers::ArrayView<const Mg::Vector3ui>& triangles,
        const Cr::Containers::ArrayView<const Mg::Vector3d>& vertices,
        bool uniform = false);

Eigen::SparseMatrix<Mg::Double> computeStiffnessMatrix(
        const Cr::Containers::ArrayView<const Mg::Vector3ui>& triangles,
        const Cr::Containers::ArrayView<const Mg::Vector3d>& vertices);
