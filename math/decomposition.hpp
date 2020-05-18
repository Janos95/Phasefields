//
// Created by janos on 17.05.20.
//

#pragma once

#include <Magnum/Math/Algorithms/Qr.h>

#include "detail.hpp"
#include "vector.hpp"
#include "matrix.hpp"

#include <Magnum/Magnum.h>
#include <Corrade/Containers/Array.h>
#include <Corrade/Containers/StridedArrayView.h>
#include <Magnum/Math/Math.h>
#include <Magnum/Math/TypeTraits.h>
#include <Magnum/Math/Functions.h>
#include <Corrade/Utility/Assert.h>
#include <Corrade/Utility/Algorithms.h>

namespace Mg = Magnum;
namespace Cr = Corrade;

template<class Scalar>
bool isZero(Scalar x) {
    return Mg::Math::abs(x) < Mg::Math::TypeTraits<Scalar>::epsilon();
}

template<class Scalar>
void house(Vector<Scalar> const& x, Vector<Scalar>& v, Scalar& beta){
    auto m = x.length();
    auto sigma = x(1, m).dot(x(1, m));
    v = Vector<Scalar>(Cr::Containers::NoInit, m);
    v(0) = Scalar{1};
    v(1, m) = x;
    if(isZero(sigma) && x(0) >= Scalar{0})
        beta = 0.;
    else if(isZero(sigma) && x(0) < Scalar{0})
        beta = -2.;
    else {
        auto mu = Mg::Math::sqrt(x(0) * x(0) + sigma);
        if(x(0) <= 0)
            v(0) = x(0) - mu;
        else
            v(0) = -sigma / (x(0) + mu);
        beta = 2 * v(0) * v(0) / (sigma + v(0) * v(0));
        v = v / v(1);
    }
    return std::make_pair(v, beta);
}


template<class Scalar>
void givens(Scalar a, Scalar b, Scalar& c, Scalar& s){
    if(isZero(b)) {
        c = Scalar{1.};
        s = Scalar{0.};
    } else {
        if(Mg::Math::abs(b) > Mg::Math::abs(a)){
           auto tao = -a/b;
           s = Mg::Math::sqrtInverted(Scalar{1.} + tao * tao);
           c = s * tao;
        } else {
            auto tao = -b/a;
            c = Mg::Math::sqrtInverted(Scalar{1.} + tao * tao);
            s = c * tao;
        }
    }
}

template<class Scalar>
struct HouseholderQR {

    HouseholderQR(Matrix<Scalar> const& A_):
       A(A_)
    {
        Mg::Int m = A.rows(), n = A.cols();
        Vector<Scalar> v;
        Scalar beta;
        for (int j = 0; j < n; ++j) {
            house(A({j,j}, {A.cols(),j + 1}), v, beta);
            A({j,j},{m,n}) = (Matrix<Scalar>::Identity() - beta * v * v.transposed()) * A({j,j}, {m,n});
            if(j < m){
                A({j+1,j}, {m, j+1}) = v(1, m - j + 1);
            }
        }
    }

    Scalar absDeterminant(){
        Scalar det{1.};
        for (int i = 0; i < A.cols(); ++i) det *= A(i,i);
        return Mg::Math::abs(det);
    }

    Vector<Scalar> solve(Vector<Scalar> const& y){
        Vector<Scalar> yp = 
    }

    Matrix<Scalar> A;
};