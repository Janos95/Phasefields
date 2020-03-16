//
// Created by janos on 03.03.20.
//

#pragma once

#include <Eigen/Core>
#include <folly/Function.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Magnum.h>
#include <string>

void writeMesh(
        const std::string& path,
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        const Eigen::VectorXd& U,
        folly::FunctionRef<Magnum::Color4(double)>);
