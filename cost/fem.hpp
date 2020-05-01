//
// Created by janos on 30.04.20.
//

#pragma once


#include <Eigen/SparseCore>

#include <Corrade/Containers/ArrayView.h>
#include <Magnum/Magnum.h>
#include <Magnum/Math/Vector3.h>

namespace Cr = Corrade;
namespace Mg = Magnum;


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
