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
        case 0: interp = (1. - w1 - w2);
            break;
        case 1: interp = w1;
            break;
        case 2: interp = w2;
            break;
        default:
            assert(false);
#ifdef NDEBUG
            __builtin_unreachable();
#endif
    }
    return interp;
}

template<typename T>
struct DoubleWell
{
    T v1, v2, v3;

    auto operator()(const T w1, const T w2) const
    {
        auto vInterp =  w1*v2 + w2 * v3 + (1. - w1 - w2) * v1;
        return 9./16. * std::pow(vInterp * vInterp - 1., 2);
    }
};

template<typename T>
struct DoubleWellGrad
{
    T v1, v2, v3;
    int i;

    auto operator()(const T w1, const T w2) const
    {
        auto vInterp =  w1*v2 + w2 * v3 + (1. - w1 - w2) * v1;
        auto pInterp = interpolateBasisFunction(w1, w2, i);

        return 9./4. * ( vInterp * vInterp - 1.) * vInterp * pInterp;
    }
};


template<class T>
struct IndicatorFunction
{
    T v1, v2, v3;

    auto operator()(const double w1, const double w2) const
    {
        double vInterp =  w1*v2 + w2 * v3 + (1. - w1 - w2) * v1;
        return 1./4. * std::pow( vInterp + 1., 2);
    }
};


template<class T>
struct IndicatorFunctionGrad
{
    T v1, v2, v3;
    int i;

    auto operator()(const T w1, const T w2) const
    {
        auto vInterp = w1*v2 + w2 * v3 + (1. - w1 - w2) * v1;
        auto pInterp = interpolateBasisFunction(w1, w2, i);

        return .5 * ( vInterp + 1.) * pInterp;
    }
};

