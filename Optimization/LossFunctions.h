#pragma once

#include "SmartEnum.h"

//class adouble;

namespace Phasefield {

SMART_ENUM(LossFunctionType, int,
           Trivial,
           Cauchy,
           Quadratic,
           Unknown)

/* type erasing wrapper around a loss function */
struct LossFunction {

    template<class T>
    LossFunction(T f);

    LossFunction(LossFunctionType::Value);

    ~LossFunction();

    LossFunction(LossFunction const&) = delete;

    LossFunction& operator=(LossFunction const&) = delete;

    LossFunction(LossFunction&& other) noexcept;

    LossFunction& operator=(LossFunction&& other) noexcept;

    void swap(LossFunction& other);

    void operator()(double const& in, double out[3]) const;

    //void operator()(adouble const& x, adouble& y) const;

    void drawSettings();

    //void (* ad)(void*, adouble const&, adouble&) = nullptr;

    void (* loss)(void*, double, double[3]) = nullptr;

    void (* draw)(void*) = nullptr;

    void (* destroy)(void*) = nullptr;

    void* erased = nullptr;

    LossFunctionType::Value lossType = LossFunctionType::Unknown;
    double weight = 1.;
};

void drawLossFunction(LossFunction&);

struct TrivialLoss {

    void operator()(double const& in, double out[3]) const;

    static LossFunctionType::Value type() { return LossFunctionType::Trivial; }

    template<class T>
    void operator()(T const& x, T& y) const {
        y = x;
    }
};

struct QuadraticLoss {

    void operator()(double const& in, double out[3]) const;

    template<class T>
    void operator()(T const& x, T& y) const {
        y = x*x;
    }

    static LossFunctionType::Value type() { return LossFunctionType::Quadratic; }
};

struct CauchyLoss {

    void operator()(double const& in, double out[3]) const;

    template<class T>
    void operator()(T const& x, T& y) const {
        y = log(1 + x);
    }

    static LossFunctionType::Value type() { return LossFunctionType::Cauchy; }
};

/* explicit template declarations */
extern template LossFunction::LossFunction(TrivialLoss);

extern template LossFunction::LossFunction(CauchyLoss);

extern template LossFunction::LossFunction(QuadraticLoss);

}