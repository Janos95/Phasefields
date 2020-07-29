//
// Created by janos on 17.05.20.
//

#pragma once

#include "matrix_expressions.hpp"
#include "tensor.hpp"

template<class Scalar>
class Matrix : public Tensor<2, Scalar> {
public:


    static IdentityExpression<Scalar> Identity(Mg::Int n) { return IdentityExpression<Scalar>{n}; }

    [[nodiscard]] constexpr inline Mg::Int rows() const noexcept { return this->m_size[0]; }

    [[nodiscard]] constexpr inline Mg::Int cols() const noexcept { return this->m_size[1]; }

    Expression<2, Scalar> inverted() {

    }

    Expression<2, Scalar> transpose() {

    }

    Expression<2, Scalar> determinant() noexcept {

    }
};

