//
// Created by janos on 06.06.20.
//

#pragma once

#include <Corrade/Utility/Macros.h>

/* In c++20 the stl <algorithm> header transitivly includes the <range> header which is huge */
#ifdef CORRADE_TARGET_LIBSTDCXX
#include <bits/stl_algobase.h>
#include <bits/stl_algo.h>
#else
#include <algorithm>
#endif
