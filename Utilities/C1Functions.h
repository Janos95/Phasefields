//
// Created by janos on 28.04.20.
//

#pragma once

#include "Types.h"

#include <Corrade/Utility/StlMath.h>
#include <Magnum/Math/Functions.h>

struct F {

    template<class T>
    T eval(const T x) const {
        if(x < a) return pow(x - a, 2);
        if(x < b) return T{0};
        return pow(x - b, 2);
    }

    template<class T>
    T grad(const T x) const {
        if(x < a) return 2*(x - a);
        if(x < b) return T{0};
        return 2*(x - b);
    }

    double a, b;
};

class W {
public:

    W(const double a, const double b) : m_a(a), m_b(b) {
        //double steepestPoint = 1./6.*(Math::sqrt(3) * Math::abs(b-a) + 3*(a+b));
        //m_scale = .01/Math::abs(grad(steepestPoint));
        m_scale = -30./pow(a - b, 5);
    }

    template<class T>
    T eval(const T& x) const {
        if(x < m_a)
            return 0.;
        if(x < m_b)
            return pow(x - m_a, 2)*pow(x - m_b, 2)*m_scale;

        return 0.;
    }

    template<class Scalar>
    inline Scalar grad(Scalar x) const {
        if(x < m_a)
            return 0.;
        if(x < m_b)
            return 2.*m_scale*((x - m_a)*pow(x - m_b, 2) + pow(x - m_a, 2)*(x - m_b));
        return 0.;
    }

private:
    const double m_a, m_b;
    double m_scale = 1.;
};

struct SmootherStep {
    template<class T>
    constexpr T eval(T x) const {
        if(x <= -1) return T{0};
        if(x <= 1.){
            x = .5*(x + 1.);
            return x*x*(3. - 2.*x);
        }
        return T{1};
    }

    template<class T>
    constexpr T grad(T x) const {
        if(x <= -1) return T{0};
        if(x <= 1.){
            x = .5*(x + 1.);
            return 3.*x*(1. - x);
        }
        return T{0};
    }
};

struct LinearChi {
    template<class T>
    constexpr T eval(T x) const {
        return 0.5*(x+1);
    }

    template<class T>
    constexpr T grad(T x) const {
        return 0.5;
    }
};

struct QuadraticChi {
    template<class T>
    constexpr T eval(T x) const {
        return 0.25*(x+1)*(x+1);
    }

    template<class T>
    constexpr T grad(T x) const {
        return 0.5*(x + 1);
    }
};

struct DoubleWell {
    template<typename T>
    auto eval(const T x) const {
        return 9./16.*(x*x - 1.)*(x*x - 1);
    }

    template<typename T>
    auto grad(const T x) const {
        return 9./4.*(x*x - 1)*x;
    }

};

struct SmoothIndicatorFunction {
    template<class T>
    auto eval(const T x) const {
        return 0.25*(x + 1.)*(x + 1.);
    }

    template<class T>
    auto grad(const T x) const {
        return .5*(x + 1.);
    }
};

struct WeightExitPenalty {

    double p = 10;

    template<class T>
    T eval(const T x) const {
        if(x < 1.)
            return p*(x - 1)*(x - 1);
        return T(0.);
    }

    template<class T>
    T grad(const T x) const {
        if(x < 1.)
            return 2*p*(x - 1);
        return T(0.);
    }
};

