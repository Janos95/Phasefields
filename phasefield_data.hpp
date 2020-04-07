//
// Created by janos on 02.04.20.
//

#pragma once

#include "drawable_data.hpp"
#include "drawables.hpp"

#include <Corrade/Containers/Array.h>

#include <Magnum/Magnum.h>
#include <Magnum/Math/Vector3.h>

#include <mutex>

struct PhasefieldData : ColorMapDrawableData {

    enum class Status: Magnum::UnsignedShort {
        NothingChanged,
        NewMesh,
        Subdivided,
        PhasefieldUpdated
    };

    Corrade::Containers::Array<Magnum::Vector3> V;
    Corrade::Containers::Array<Magnum::UnsignedInt> F;
    Corrade::Containers::Array<double> phasefield;
    Magnum::Trade::MeshData original{Magnum::MeshPrimitive::Points, 0};
    Status status = Status::NothingChanged;
    Drawable* drawable = nullptr;
    DrawableType type;
    std::mutex mutex;
};