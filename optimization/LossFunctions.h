#pragma once

#include "SmartEnum.h"

class adouble;

SMART_ENUM(LossFunctionType, int,
           TrivialLoss,
           ScaledLoss,
           CauchyLoss,
           QuadraticLoss,
           ComposedLoss)

/* type erasing wrapper around a loss function */
struct LossFunction {

    template<class T>
    LossFunction(T f);

    explicit LossFunction(LossFunctionType::Value);

    ~LossFunction();

    LossFunction(LossFunction const&) = delete;

    LossFunction& operator=(LossFunction const&) = delete;

    LossFunction(LossFunction&& other) noexcept;

    LossFunction& operator=(LossFunction&& other) noexcept;

    void swap(LossFunction& other);

    friend void swap(LossFunction& f1, LossFunction& f2);

    void operator()(double const& in, double out[3]) const;

    void operator()(adouble const& x, adouble& y) const;

    void drawSettings();

    void (* ad)(void*, adouble const&, adouble&) = nullptr;

    void (* loss)(void*, double, double[3]) = nullptr;

    void (* draw)(void*) = nullptr;

    void (* destroy)(void*) = nullptr;

    void* erased = nullptr;

    LossFunctionType::Value lossType;
};

void drawLossFunction(LossFunction&);

struct TrivialLoss {

    void operator()(double const& in, double out[3]) const;

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
};

struct CauchyLoss {

    void operator()(double const& in, double out[3]) const;

    template<class T>
    void operator()(T const& x, T& y) const {
        y = log(1 + x);
    }

};

struct ScaledLoss {
    explicit ScaledLoss(LossFunction f_, double s_);

    void operator()(double const& in, double out[3]) const;

    template<class T>
    void operator()(T const& x, T& y) const {
        f(x, y);
        y *= s;
    }

    LossFunction f;
    double s;
};

struct ComposedLoss {
    explicit ComposedLoss(LossFunction g_, LossFunction f_);

    void operator()(double const& in, double out[3]) const;

    template<class T>
    void operator()(T const& x, T& y) const {
        T temp;
        f(x, temp);
        g(temp, y);
    }

    LossFunction g, f;
};

/* explicit template declarations */
extern template LossFunction::LossFunction(TrivialLoss);

extern template LossFunction::LossFunction(ScaledLoss);

extern template LossFunction::LossFunction(CauchyLoss);

extern template LossFunction::LossFunction(QuadraticLoss);

extern template LossFunction::LossFunction(ComposedLoss);
