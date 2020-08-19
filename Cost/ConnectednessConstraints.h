//
// Created by janos on 21.05.20.
//

#pragma once

#include "small_array.hpp"
#include "Paths.h"
#include "Functional.h"
#include "Dijkstra.hpp"
#include "StoppingCriteria.h"

#include <Corrade/Containers/Array.h>

namespace Mn = Magnum;
namespace Cr = Corrade;

struct ConnectednessConstraint;
template<>
Functional::Functional<ConnectednessConstraint>;
