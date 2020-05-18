//
// Created by janos on 17.05.20.
//

#pragma once

template<unsigned dimensions, class Scalar>
struct Expression {

    using SizeType = typename Cr::Containers::StridedArrayView<dimensions, Scalar>::Size;

    virtual Tensor<dimensions, Scalar> evaluate() = 0;
    virtual Tensor<dimensions, Scalar> optimize() {};
    virtual Cr::Containers::ArrayView<Cr::Containers::Pointer<Expression>> children() = 0;
    virtual std::size_t cost() = 0;
    virtual SizeType size() = 0;
};

template<class Scalar>
struct DimErasedExpression {
    unsigned dimension;
};
