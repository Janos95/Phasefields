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

struct ParametricSmoothStep {
    double edge0, edge1;

    template<class Scalar>
    Scalar clamp(Scalar x) {
        if (x < 0)
            x = 0;
        if (x > 1)
            x = 1;
        return x;
    }

    template<class Scalar>
    Scalar eval(Scalar x) {
        Scalar xScaled = (x - edge0) / (edge1 - edge0);
        x = clamp(xScaled);
        return x * x * (3 - 2 * x);
    }

    template<class Scalar>
    Scalar grad(Scalar x) {
        Scalar scale = 1./(edge1 - edge0);
        Scalar xScaled = (x - edge0)*scale;
        x = clamp(xScaled);
        return -6*(-1 + x)*x*scale;
    }
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

struct QuadraticChiMirrored {
    template<class T>
    constexpr T eval(T x) const {
        return 0.25*(x-1)*(x-1);
    }

    template<class T>
    constexpr T grad(T x) const {
        return 0.5*(x - 1);
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

    const double eps = 1e-6;

    template<class T>
    T eval(const T x) const {
        return 1./(x + eps);
    }

    template<class T>
    T grad(const T x) const {
        return -1./Math::pow<2>(x+eps);
    }
};

struct Quadratic {

    template<class T>
    constexpr T eval(const T x) const {
        return x*x;
    }

    template<class T>
    constexpr T grad(const T x) const {
        return 2*x;
    }
};

struct Quartic {

    template<class T>
    constexpr T eval(const T x) const {
        if(x < -1) return 1;
        if(x < 1) return x*x*(2-x*x);
        return 1;
    }

    template<class T>
    constexpr T grad(const T x) const {
        if(x < -1) return 0;
        if(x < 1) return 2*x*(1-2*x*x);
        return 0;
    }
};

