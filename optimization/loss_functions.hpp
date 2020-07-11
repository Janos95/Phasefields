#pragma once

class adouble;

/* type erasing wrapper around a loss function */
struct LossFunction {

    template<class T>
    LossFunction(T f);

    ~LossFunction();

    LossFunction(LossFunction const&) = delete;

    LossFunction& operator=(LossFunction const&) = delete;

    LossFunction(LossFunction&& other) noexcept;

    LossFunction& operator=(LossFunction&& other) noexcept;

    friend void swap(LossFunction& f1, LossFunction& f2);

    void operator()(double const& in, double out[3]);

    void operator()(adouble const& x, adouble& y);

    void (*ad)(void*, adouble const&, adouble&) = nullptr;
    void (*loss)(void*, double, double[3]) = nullptr;
    void (*destroy)(void*) = nullptr;
    void* erased = nullptr;
};

struct TrivialLoss {

    void operator()(double const& in, double out[3]);

    template<class T>
    void operator()(T const& x, T& y){
        y = x;
    }
};

struct QuadraticLoss {

    void operator()(double const& in, double out[3]);

    template<class T>
    void operator()(T const& x, T& y){
        y = x * x;
    }
};

struct CauchyLoss  {

    void operator()(double const& in, double out[3]);

    template<class T>
    void operator()(T const& x, T& y){
        y = log(1 + x);
    }

};

struct ScaledLoss  {
    explicit ScaledLoss(LossFunction f_, double s_);

    void operator()(double const& in, double out[3]);

    template<class T>
    void operator()(T const& x, T& y){
        y = f(x) * s;
    }

    LossFunction f;
    double s;
};

struct ComposedLoss  {
    explicit ComposedLoss(LossFunction g_, LossFunction f_);

    void operator()(double const& in, double out[3]);

    template<class T>
    T operator()(T const& x, T& y){
        y = g(f(x));
    }

    LossFunction g, f;
};

/* explicit template declarations */
extern template LossFunction::LossFunction(TrivialLoss);
extern template LossFunction::LossFunction(ScaledLoss);
extern template LossFunction::LossFunction(CauchyLoss);
extern template LossFunction::LossFunction(QuadraticLoss);
extern template LossFunction::LossFunction(ComposedLoss);
