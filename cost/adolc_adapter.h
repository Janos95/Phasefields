//
// Created by janos on 09.06.20.
//

#pragma once

#include "functional.hpp"
#include "types.hpp"

#include <Corrade/Containers/Pointer.h>
#include <adolc/adouble.h>

static int g_numTags = 0;

template<template <class> class F>
struct AdolcAdapter final : FunctionalD {

    template<class... Args>
    explicit AdolcAdapter(Args const&...);
    ~AdolcAdapter() override;

    Pointer<F<adouble>> ad;
    Pointer<F<double>> plain;

    int tag;

    bool evaluate(double const* parameter, double* cost, double** jacobian) const override;
    [[nodiscard]] uint32_t numParameters() const override;
    [[nodiscard]] uint32_t numResiduals() const override;
};
