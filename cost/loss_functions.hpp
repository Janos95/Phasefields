#pragma once

#include <Corrade/Containers/Pointer.h>
#include <Magnum/Magnum.h>

namespace Cr = Corrade;
namespace Mg = Magnum;

enum class LossType : Mg::UnsignedInt {
    Trivial = 0,
    Huber = 1,
    SoftOne = 2,
    Cauchy = 3,
    Arctan = 4,
    Tolerant = 5,
    Tukey = 6,
    Composed = 7,
    Scaled = 8
};

struct LossFunction {
  virtual ~LossFunction() = default;

  // For a residual vector with squared 2-norm 'sq_norm', this method
  // is required to fill in the value and derivatives of the loss
  // function (rho in this example):
  //
  //   out[0] = rho(sq_norm),
  //   out[1] = rho'(sq_norm),
  //   out[2] = rho''(sq_norm),
  //
  // Here the convention is that the contribution of a term to the
  // cost function is given by 1/2 rho(s),  where
  //
  //   s = ||residuals||^2.
  //
  // Calling the method with a negative value of 's' is an error and
  // the implementations are not required to handle that case.
  //
  // Most sane choices of rho() satisfy:
  //
  //   rho(0) = 0,
  //   rho'(0) = 1,
  //   rho'(s) < 1 in outlier region,
  //   rho''(s) < 0 in outlier region,
  //
  // so that they mimic the least squares cost for small residuals.
  virtual void Evaluate(double sq_norm, double out[3]) const = 0;

  LossType type;
};

// Some common implementations follow below.
//
// Note: in the region of interest (i.e. s < 3) we have:
//   TrivialLoss >= HuberLoss >= SoftLOneLoss >= CauchyLoss

// This corresponds to no robustification.
//
//   rho(s) = s
//
// At s = 0: rho = [0, 1, 0].
//
// It is not normally necessary to use this, as passing NULL for the
// loss function when building the problem accomplishes the same
// thing.
struct TrivialLoss final : LossFunction {
    TrivialLoss()  { type = LossType::Trivial; };
  void Evaluate(double, double*) const override;
};

// Scaling
// -------
// Given one robustifier
//   s -> rho(s)
// one can change the length scale at which robustification takes
// place, by adding a scale factor 'a' as follows:
//
//   s -> a^2 rho(s / a^2).
//
// The first and second derivatives are:
//
//   s -> rho'(s / a^2),
//   s -> (1 / a^2) rho''(s / a^2),
//
// but the behaviour near s = 0 is the same as the original function,
// i.e.
//
//   rho(s) = s + higher order terms,
//   a^2 rho(s / a^2) = s + higher order terms.
//
// The scalar 'a' should be positive.
//
// The reason for the appearance of squaring is that 'a' is in the
// units of the residual vector norm whereas 's' is a squared
// norm. For applications it is more convenient to specify 'a' than
// its square. The commonly used robustifiers below are described in
// un-scaled format (a = 1) but their implementations work for any
// non-zero value of 'a'.

// Huber.
//
//   rho(s) = s               for s <= 1,
//   rho(s) = 2 sqrt(s) - 1   for s >= 1.
//
// At s = 0: rho = [0, 1, 0].
//
// The scaling parameter 'a' corresponds to 'delta' on this page:
//   http://en.wikipedia.org/wiki/Huber_Loss_Function
struct HuberLoss final : LossFunction {
  explicit HuberLoss(double a_) : a(a_), b(a_ * a_) { type = LossType::Huber; }
  void Evaluate(double, double*) const override;

  const double a;
  // b = a^2.
  const double b;
};

// Soft L1, similar to Huber but smooth.
//
//   rho(s) = 2 (sqrt(1 + s) - 1).
//
// At s = 0: rho = [0, 1, -1 / (2 * a^2)].
struct SoftLOneLoss final : LossFunction {
  explicit SoftLOneLoss(double a_) : b(a_ * a_), c(1.f / b) { type = LossType::SoftOne; }
  void Evaluate(double, double*) const override;

  // b = a^2.
  const double b;
  // c = 1 / a^2.
  const double c;
};

// Inspired by the Cauchy distribution
//
//   rho(s) = log(1 + s).
//
// At s = 0: rho = [0, 1, -1 / a^2].
struct CauchyLoss final : LossFunction {
  explicit CauchyLoss(double a_) : b(a_ * a_), c(1.f / b) { type = LossType::Cauchy; }
  void Evaluate(double, double*) const override;

  // b = a^2.
  const double b;
  // c = 1 / a^2.
  const double c;
};

// Loss that is capped beyond a certain level using the arc-tangent function.
// The scaling parameter 'a' determines the level where falloff occurs.
// For costs much smaller than 'a', the loss function is linear and behaves like
// TrivialLoss, and for values much larger than 'a' the value asymptotically
// approaches the constant value of a * PI / 2.
//
//   rho(s) = a atan(s / a).
//
// At s = 0: rho = [0, 1, 0].
struct ArctanLoss final : LossFunction {
  explicit ArctanLoss(double a_) : a(a_), b(1.f / (a_ * a_)) { type = LossType::Arctan; }
  void Evaluate(double, double*) const override;

  const double a;
  // b = 1 / a^2.
  const double b;
};

// Loss function that maps to approximately zero cost in a range around the
// origin, and reverts to linear in error (quadratic in cost) beyond this range.
// The tolerance parameter 'a' sets the nominal point at which the
// transition occurs, and the transition size parameter 'b' sets the nominal
// distance over which most of the transition occurs. Both a and b must be
// greater than zero, and typically b will be set to a fraction of a.
// The slope rho'[s] varies smoothly from about 0 at s <= a - b to
// about 1 at s >= a + b.
//
// The term is computed as:
//
//   rho(s) = b log(1 + exp((s - a) / b)) - c0.
//
// where c0 is chosen so that rho(0) == 0
//
//   c0 = b log(1 + exp(-a / b)
//
// This has the following useful properties:
//
//   rho(s) == 0               for s = 0
//   rho'(s) ~= 0              for s << a - b
//   rho'(s) ~= 1              for s >> a + b
//   rho''(s) > 0              for all s
//
// In addition, all derivatives are continuous, and the curvature is
// concentrated in the range a - b to a + b.
//
// At s = 0: rho = [0, ~0, ~0].
struct TolerantLoss final : LossFunction {
  explicit TolerantLoss(double a_, double b_);
  void Evaluate(double, double*) const override;

  const double a, b, c;
};

// This is the Tukey biweight loss function which aggressively
// attempts to suppress large errors.
//
// The term is computed as follows where the equations are scaled by a
// factor of 2 because the cost function is given by 1/2 rho(s):
//
//   rho(s) = a^2 / 3 * (1 - (1 - s / a^2)^3 )   for s <= a^2,
//   rho(s) = a^2 / 3                          for s > a^2.
//
// At s = 0: rho = [0, 1, -2 / a^2]
struct TukeyLoss : LossFunction {
  explicit TukeyLoss(double a) : a_squared(a * a) { type = LossType::Tukey; }
  void Evaluate(double, double*) const override;

  const double a_squared;
};

// Composition of two loss functions.  The error is the result of first
// evaluating g followed by f to yield the composition f(g(s)).
// The loss functions must not be NULL.
struct ComposedLoss final : LossFunction {
  explicit ComposedLoss(Cr::Containers::Pointer<LossFunction> f, Cr::Containers::Pointer<LossFunction> g);

  void Evaluate(double, double*) const override;

  Cr::Containers::Pointer<LossFunction> f, g;
};

// The discussion above has to do with length scaling: it affects the space
// in which s is measured. Sometimes you want to simply scale the output
// value of the robustifier. For example, you might want to weight
// different error terms differently (e.g., weight pixel reprojection
// errors differently from terrain errors).
//
// If rho is the wrapped robustifier, then this simply outputs
// s -> a * rho(s)
//
// The first and second derivatives are, not surprisingly
// s -> a * rho'(s)
// s -> a * rho''(s)
//
// Since we treat the a NULL Loss function as the Identity loss
// function, rho = NULL is a valid input and will result in the input
// being scaled by a. This provides a simple way of implementing a
// scaled ResidualBlock.
struct ScaledLoss : LossFunction {

  ScaledLoss(Cr::Containers::Pointer<LossFunction> f_, double a_)
      : f(std::move(f_)), a(a_) { type = LossType::Scaled; }

  void Evaluate(double, double*) const override;

  Cr::Containers::Pointer<LossFunction> f;
  const double a;
};



