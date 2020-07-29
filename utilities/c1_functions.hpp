//
// Created by janos on 28.04.20.
//

#pragma once

#include <Corrade/Utility/StlMath.h>
#include <Magnum/Math/Functions.h>

template<class T>
struct F {
    T eval(const T x) const {
        if(x < a) return pow(x - a, 2);
        if(x < b) return T{0};
        return pow(x - b, 2);
    }

    T grad(const T x) const {
        if(x < a) return 2*(x - a);
        if(x < b) return T{0};
        return 2*(x - b);
    }

    T a, b;
};

class W {
public:

    W(const double a, const double b) : m_a(a), m_b(b), m_c3(-30./std::pow(a - b, 5)) {
    }

    template<class T>
    T operator()(const T& x) const {
        if(x < m_a)
            return 0.;
        if(x < m_b)
            return pow(x - m_a, 2)*pow(x - m_b, 2)*m_c3;

        return 0.;
    }

private:
    const double m_a, m_b, m_c3;
};

class WGrad {
public:

    WGrad(const double a, const double b) : m_a(a), m_b(b), m_c3(-30./std::pow(a - b, 5)) {
    }

    inline double operator()(const double x) const {
        if(x < m_a)
            return 0.;
        if(x < m_b)
            return 2.*((x - m_a)*pow(x - m_b, 2) + pow(x - m_a, 2)*(x - m_b))*m_c3;
        return 0.;
    }

private:
    const double m_a, m_b, m_c3;
};

struct DoubleWell {
    template<typename T>
    auto eval(const T x) const {
        return 9./16.*std::pow(x*x - 1., 2);
    }

    template<typename T>
    auto grad(const T x) const {
        return 9./4.*(x*x - 1)*x;
    }

};


struct SmootherStep {

    template<class T>
    T eval(T x) const {
        if(x <= -1) return T{0};
        if(x <= 1.){
            x = .5*(x + 1.);
            return x*x*(3. - 2.*x);
        }
        return T{1};
    }

    template<class T>
    T grad(T x) const {
        if(x <= -1) return T{0};
        if(x <= 1.){
            x = .5*(x + 1.);
            return 3.*x*(1. - x);
        }
        return T{0};
    }
};


template<class T>
struct Indicator {
    auto eval(const T x) const {
        return 0.25*std::pow(x + 1., 2);
    }

    auto grad(const T x) const {
        return .5*(x + 1.);
    }
};

