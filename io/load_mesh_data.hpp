//
// Created by janos on 2/8/20.
//

#pragma once

#include <Magnum/Trade/MeshData.h>

#include <Eigen/Core>

#include <string>

Magnum::Trade::MeshData loadMeshData(const std::string& path);

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi> toEigen(Magnum::Trade::MeshData const& mesh);