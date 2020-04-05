//
// Created by janos on 02.04.20.
//

#pragma once

#include "drawable_data.hpp"

#include <Corrade/Containers/Array.h>

#include <Magnum/Magnum.h>
#include <Magnum/Math/Vector3.h>

#include <mutex>

struct PhasefieldData : DrawableData {

    enum class Status: Magnum::UnsignedShort {
        NewMesh,
        Subdivided,
        PhasefieldUpdated
    };

    Corrade::Containers::Array<Magnum::Vector3d> vertices;
    Corrade::Containers::Array<Int> faces;
    Corrade::Containers::Array<double> phasefield;
    Magnum::Trade::MeshData original;
    Status status;
    std::mutex mutex;
};