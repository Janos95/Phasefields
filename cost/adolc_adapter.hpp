//
// Created by janos on 09.06.20.
//

#include "adolc_adapter.h"

#include <adolc/adolc.h>

template<template <class> class F>
template<class... Args>
AdolcAdapter<F>::AdolcAdapter(Args const&... args) :
    ad(args...), plain(args...), tag(g_numTags++)
{
}

template<template <class> class F>
AdolcAdapter<F>::~AdolcAdapter(){ --g_numTags; }

template<template <class> class F>
bool AdolcAdapter<F>::evaluate(double const* parameter, double* cost, double** jac) const {
    if(jac){
        trace_on(tag);
        {
            Array<adouble> xs(numParameters());
            Array<adouble> ys(numResiduals());
            for (uint32_t i = 0; i < numParameters(); ++i)
                xs[i] <<= parameter[i];
            ad->evaluate(xs.data(), ys.data(), nullptr);
            for (uint32_t i = 0; i < ys.size(); ++i)
                ys[i] >>= cost[i];
        }
        trace_off();
        jacobian(tag, numResiduals(), numParameters(), parameter, jac);
    } else {
        plain->evaluate(parameter, cost, nullptr);
    }
    return true;
}

template<template <class> class F>
uint32_t AdolcAdapter<F>::numParameters() const {
    return plain->numParameters();
}

template<template <class> class F>
uint32_t AdolcAdapter<F>::numResiduals() const {
    return ad->numResiduals();
}
