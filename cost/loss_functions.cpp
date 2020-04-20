//
// Created by janos on 20.04.20.
// from ceres/loss_function.h

#include "loss_functions.hpp"

#include <cmath>

using namespace Corrade;

void TrivialLoss::Evaluate(double s, double rho[3]) const {
    rho[0] = s;
    rho[1] = 1.0;
    rho[2] = 0.0;
}

void HuberLoss::Evaluate(double s, double rho[3]) const {
    if (s > b) {
        // Outlier region.
        // 'r' is always positive.
        const double r = sqrt(s);
        rho[0] = 2.0 * a * r - b;
        rho[1] = std::max(std::numeric_limits<double>::min(), a / r);
        rho[2] = - rho[1] / (2.0 * s);
    } else {
        // Inlier region.
        rho[0] = s;
        rho[1] = 1.0;
        rho[2] = 0.0;
    }
}

void SoftLOneLoss::Evaluate(double s, double rho[3]) const {
    const double sum = 1.0 + s * c;
    const double tmp = sqrt(sum);
    // 'sum' and 'tmp' are always positive, assuming that 's' is.
    rho[0] = 2.0 * b * (tmp - 1.0);
    rho[1] = std::max(std::numeric_limits<double>::min(), 1.0 / tmp);
    rho[2] = - (c * rho[1]) / (2.0 * sum);
}

void CauchyLoss::Evaluate(double s, double rho[3]) const {
    const double sum = 1.0 + s * c;
    const double inv = 1.0 / sum;
    // 'sum' and 'inv' are always positive, assuming that 's' is.
    rho[0] = b * log(sum);
    rho[1] = std::max(std::numeric_limits<double>::min(), inv);
    rho[2] = - c * (inv * inv);
}

void ArctanLoss::Evaluate(double s, double rho[3]) const {
    const double sum = 1 + s * s * b;
    const double inv = 1 / sum;
    // 'sum' and 'inv' are always positive.
    rho[0] = a * atan2(s, a);
    rho[1] = std::max(std::numeric_limits<double>::min(), inv);
    rho[2] = -2.0 * s * b * (inv * inv);
}

TolerantLoss::TolerantLoss(double a_, double b_)
        : a(a_),
          b(b_),
          c(b_ * log(1.0 + exp(-a_ / b_))) {
    type = LossType::Tolerant;
    CORRADE_ASSERT(a >= 0.0, "tolerant loss: a needs to be greater or equal 0",);
    CORRADE_ASSERT(b > 0.0, "tolerant loss: b needs to bigger than 0", );
}

void TolerantLoss::Evaluate(double s, double rho[3]) const {
    const double x = (s - a) / b;
    // The basic equation is rho[0] = b ln(1 + e^x).  However, if e^x is too
    // large, it will overflow.  Since numerically 1 + e^x == e^x when the
    // x is greater than about ln(2^53) for doubles, beyond this threshold
    // we substitute x for ln(1 + e^x) as a numerically equivalent approximation.
    static const double kLog2Pow53 = 36.7;  // ln(MathLimits<double>::kEpsilon).
    if (x > kLog2Pow53) {
        rho[0] = s - a - c;
        rho[1] = 1.0;
        rho[2] = 0.0;
    } else {
        const double e_x = exp(x);
        rho[0] = b * log(1.0 + e_x) - c;
        rho[1] = std::max(std::numeric_limits<double>::min(), e_x / (1.0 + e_x));
        rho[2] = 0.5 / (b * (1.0 + cosh(x)));
    }
}

void TukeyLoss::Evaluate(double s, double* rho) const {
    if (s <= a_squared) {
        // Inlier region.
        const double value = 1.0 - s / a_squared;
        const double value_sq = value * value;
        rho[0] = a_squared / 3.0 * (1.0 - value_sq * value);
        rho[1] = value_sq;
        rho[2] = -2.0 / a_squared * value;
    } else {
        // Outlier region.
        rho[0] = a_squared / 3.0;
        rho[1] = 0.0;
        rho[2] = 0.0;
    }
}

ComposedLoss::ComposedLoss(Containers::Pointer<LossFunction> f_, Containers::Pointer<LossFunction> g_)
        : f(std::move(f_)),
          g(std::move(g_)){
    type = LossType::Composed;
    CORRADE_ASSERT(f, "Composed Loss: f must not be nullptr",);
    CORRADE_ASSERT(g, "Composed Loss: g must not be nullptr",);
}


void ComposedLoss::Evaluate(double s, double rho[3]) const {
    double rho_f[3], rho_g[3];
    g->Evaluate(s, rho_g);
    f->Evaluate(rho_g[0], rho_f);
    rho[0] = rho_f[0];
    // f'(g(s)) * g'(s).
    rho[1] = rho_f[1] * rho_g[1];
    // f''(g(s)) * g'(s) * g'(s) + f'(g(s)) * g''(s).
    rho[2] = rho_f[2] * rho_g[1] * rho_g[1] + rho_f[1] * rho_g[2];
}

void ScaledLoss::Evaluate(double s, double rho[3]) const {
    if (!f) {
        rho[0] = a * s;
        rho[1] = a;
        rho[2] = 0.0;
    } else {
        f->Evaluate(s, rho);
        rho[0] *= a;
        rho[1] *= a;
        rho[2] *= a;
    }
}

