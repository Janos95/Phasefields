//
// Created by janos on 20.04.20.
//

#pragma once

#include "loss_functions.hpp"
#include "shared_ressource.hpp"

#include <Corrade/Containers/Pointer.h>
#include <string>

namespace Cr = Corrade;
namespace Mg = Magnum;

enum class FunctionalType : Mg::UnsignedInt {
    Undefined = 0,
    DirichletEnergy = 1,
    DoubleWellPotential = 2,
    Area1 = 3,
    Connectedness = 4,
    Area2 = 5,
};



struct Functional{
    std::string name;
    Cr::Containers::Pointer<LossFunction> loss = Cr::Containers::pointer<TrivialLoss>();
    SharedRessource<Mg::Double> scaling = nullptr;
    FunctionalType type;

    virtual bool evaluate(double const* parameter, double* cost, double* jacobian) const = 0;
    virtual int numParameters() const = 0;
    virtual ~Functional() = default;
};
