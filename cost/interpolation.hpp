//
// Created by janos on 27.11.19.
//

#pragma once

#include <cmath>
#include <cassert>

template<typename T>
auto interpolateBasisFunction(const T w1, const T w2, const int basis)
{
    T interp;
    switch(basis)
    {
        case 0: return (T(1.) - w1 - w2);
        case 1: return w1;
        case 2: return w2;
        default:
            assert(false);
#ifdef NDEBUG
            __builtin_unreachable();
#endif
    }
}


