//
// Created by janos on 03.03.20.
//

#pragma once

#include <Eigen/Core>

#include <Corrade/Containers/ArrayView.h>
#include <Magnum/Math/Vector3.h>
#include <Magnum/Magnum.h>


void writeMesh(
        const std::string&,
        Corrade::Containers::ArrayView<Magnum::Vector3> varr,
        Corrade::Containers::ArrayView<Magnum::UnsignedInt> farr,
        const Eigen::VectorXd&,
        bool = false);