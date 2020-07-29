//
// Created by janos on 21.05.20.
//

#pragma once

#include "small_array.hpp"
#include "paths.hpp"
#include "Functional.h"
#include "dijkstra.hpp"
#include "stopping_criteria.hpp"

#include <Corrade/Containers/Array.h>

namespace Mn = Magnum;
namespace Cr = Corrade;

struct ConnectednessConstraint;
template<>
Functional::Functional<ConnectednessConstraint>;
