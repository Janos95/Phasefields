//
// Created by janos on 17.05.20.
//

#include "detail.hpp"
#include "tensor.hpp"
#include "abstract_expression.hpp"

#include <Magnum/Magnum.h>
#include <Corrade/Containers/Array.h>
#include <Corrade/Containers/StridedArrayView.h>
#include <Magnum/Math/Math.h>
#include <Magnum/Math/TypeTraits.h>
#include <Magnum/Math/Functions.h>
#include <Corrade/Utility/Assert.h>
#include <Corrade/Utility/Algorithms.h>

template<unsigned dimensions, class Scalar>
struct ScalarExpression : Expression<dimensions, Scalar> {

    std::size_t cost(){
        return
    }

    Scalar alpha;

};

template<unsigned dimensions, class Scalar>
struct ReductionExpression : Expression<dimensions, Scalar> {
    using Size = Cr::Containers::StridedDimensions<dimensions + 1, std::size_t>;

    Tensor<dimensions, Scalar> evaluate() override {
        Tensor<dimensions, Scalar> T;
        Tensor<dimensions + 1, Scalar> S = expr->evaluate();

        for (int i = 0; i < dim; ++i) {
            Size lower{}, upper = S.size();
            lower[i]






            T(lower, upper) =
        }
        return T;
    }

    std::size_t cost(){


    }

    Corrade::Containers::Pointer<Expression<dimensions + 1, Scalar>> expr;
    int dim; //dimension along which to reduce
};
