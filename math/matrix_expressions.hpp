//
// Created by janos on 17.05.20.
//

#pragma once

#include "detail.hpp"
#include "abstract_expression.hpp"
#include "matrix.hpp"
#include "sparse_matrix.hpp"

#include <Magnum/Magnum.h>
#include <Corrade/Containers/Array.h>
#include <Corrade/Containers/StridedArrayView.h>
#include <Magnum/Math/Math.h>
#include <Magnum/Math/TypeTraits.h>
#include <Magnum/Math/Functions.h>
#include <Corrade/Utility/Assert.h>
#include <Corrade/Utility/Algorithms.h>

#include <cblas.h>

namespace Mg = Magnum;
namespace Cr = Corrade;

template<class Scalar>
struct GemmExpression : Expression<2, Scalar> {

    Matrix<Scalar> evaluate() override {
        Matrix<Scalar> A,B,C;
        A = exprA->evaluate();
        B = exprB->evaluate();
        C = exprC->evaluate();
        auto rows = A.rows();
        auto n = A.cols();
        auto cols = B.rows();
        CORRADE_ASSERT(n == B.rows(), "Matrices not compatible for product", {});

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                Scalar p{0.};
                for (int k = 0; k < n; ++k) {
                    p += A(i, k) * B(k, j);
                }
                C(i,j) = beta * C(i,j) + p;
            }
        }
        return C;
    }

    std::size_t cost(){
        auto sA = exprA->size();
        auto sB = exprB->size();
        return sA[0] * sA[1] * sB[1];
    }

    Corrade::Containers::Pointer<Expression<2, Scalar>> exprA, exprB, exprC;
    Scalar alpha, beta;
};


template<class Scalar>
struct BlasGemmExpression : Expression<2, Scalar> {

    Matrix<Scalar> evaluate() override {
        Matrix<Scalar> A,B,C;
        A = exprA->evaluate();
        B = exprB->evaluate();
        C = exprC->evaluate();
        auto rows = A.rows();
        auto n = A.cols();
        auto cols = B.rows();
    }

    std::size_t cost(){
        auto sA = exprA->size();
        auto sB = exprB->size();
        return sA[0] * sA[1] * sB[1];
    }

    Corrade::Containers::Pointer<Expression<2, Scalar>> exprA, exprB, exprC;
    Scalar alpha, beta;
};

template<class Scalar>
struct IdentityExpression : Expression<2, Scalar> {

    Matrix<Scalar> evaluate() override {
        Matrix<Scalar> A(n, n);
        for (int i = 0; i < n; ++i) A(i,i) = Scalar{1.};
        return A;
    }

    std::size_t cost(){
        return n*n;
    }

    Mg::Int n;
};
