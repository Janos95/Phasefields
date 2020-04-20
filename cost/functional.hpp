//
// Created by janos on 20.04.20.
//

#pragma once

#include "loss_functions.hpp"

#include <Corrade/Containers/Pointer.h>
#include <string>

namespace Cr = Corrade;
namespace Mg = Magnum;

enum class FunctionalType : Mg::UnsignedInt {
    Undefined = 0,
    DirichletEnergy = 1,
    DoubleWellPotential = 2,
    Area = 3,
    Connectedness = 4,
};

struct Functional{
    std::string name;
    Cr::Containers::Pointer<LossFunction> loss;
    FunctionalType type;

    virtual bool evaluate(double const* parameter, double const* jacobian) = 0;
    virtual int numParameters() = 0;
};
