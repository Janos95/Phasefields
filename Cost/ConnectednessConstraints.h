//
// Created by janos on 21.05.20.
//

#pragma once

#include "SmallArray.h"
#include "Paths.h"
#include "Functional.h"
#include "Dijkstra.h"
#include "StoppingCriteria.h"

#include <Corrade/Containers/Array.h>

namespace Mn = Magnum;
namespace Cr = Corrade;

struct ConnectednessConstraint;
template<>
Functional::Functional<ConnectednessConstraint>;
