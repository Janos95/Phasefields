//
// Created by janos on 27.11.19.
//

#pragma once

#include <Eigen/Core>

auto initAlongZAxis(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::VectorXd& U);

void initRandomNormal(Eigen::VectorXd& V);

