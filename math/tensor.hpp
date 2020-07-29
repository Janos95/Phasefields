//
// Created by janos on 16.05.20.
//

#pragma once

#include "detail.hpp"
#include "abstract_expression.hpp"

#include <Magnum/Magnum.h>
#include <Corrade/Containers/Array.h>
#include <Corrade/Containers/StridedArrayView.h>
#include <Magnum/Math/Math.h>
#include <Magnum/Math/TypeTraits.h>
#include <Magnum/Math/Functions.h>
#include <Corrade/Utility/Assert.h>
#include <Corrade/Utility/Algorithms.h>

#include <cstdlib>
#include <cblas.h>

namespace Mg = Magnum;
namespace Cr = Corrade;


template<unsigned dimensions, typename Scalar>
void compile(Expression<dimensions, Scalar>& expression) {
    expression.optimize();
    for(auto const& child : expression.children())
        compile(*child);
    expression.compiled = true;
}

template<unsigned dimensions, class Scalar>
class Map {
public:

    using ArrayView = Cr::Containers::StridedArrayView<dimensions, Scalar>;

    template<unsigned dim>
    using Size = Cr::Containers::StridedDimensions<dim, Scalar>;

    explicit Map(Cr::Containers::StridedArrayView<dimensions, Scalar> const& view) :
            m_view(view) {
    }

    template<class First, class... Next>
    Scalar& operator()(First first, Next... next) {
        static_assert(sizeof...(Next) + 1 == dimensions);
        auto data = static_cast<Scalar*>(this->_data);
        return data[index_into(std::make_index_sequence<dimensions>{}, this->_stride, first, next...)];
    };

    template<unsigned newDimensions = dimensions>
    Map operator()(Size<newDimensions> const& lower, Size<newDimensions> const& upper) {
        return Map{m_view.slice(lower, upper)};
    }

private:

    ArrayView m_view;
};

template<unsigned dimensions, typename Scalar>
class Tensor {

    using ArrayView = Cr::Containers::StridedArrayView<dimensions, Scalar>;
    using Size = typename ArrayView::Size;

    template<class TagType, class First, class... Next>
    explicit Tensor(TagType, First first, Next... next) :
            m_data(TagType{}, length(Size(first, next...))),
            m_view({m_data.data(), m_data.size()}, Size(first, next...)) {
    }

    template<class First, class... Next>
    explicit Tensor(First first, Next... next) :
            m_data(Cr::Containers::ValueInit, length(Size(first, next...))),
            m_view({m_data.data(), m_data.size()}, Size(first, next...)) {
    }

    explicit Tensor(Expression<dimensions, Scalar>& expression) {
        compile(expression);
        *this = expression.evaluate();
    }

    Tensor& operator=(Tensor other) noexcept {
        other.swap(*this);
        return *this;
    }

    Tensor(Tensor const& other) :
            m_data(choose_tag_t<Scalar>{}, other.size()),
            m_view(other.m_view) {
        Cr::Utility::copy(other.m_data, m_data);
    }

    Tensor(Tensor&& other) noexcept:
            m_data(std::move(other.m_data)),
            m_view(other.m_view) {
    }

    void swap(Tensor& other) {
        auto data = std::move(other.m_data);
        auto view = static_cast<ArrayView>(*this);
        other.m_data = std::move(m_data);
        m_data = std::move(data);
    }

    Map<dimensions, Scalar> operator()(Size const& lower, Size const& upper) {
        return Map<dimensions, Scalar>{m_view.slice(lower, upper)};
    }

    [[nodiscard]] constexpr inline Size size() const noexcept { return this->size(); }

protected:
    Cr::Containers::Array<Scalar> m_data;
    Cr::Containers::StridedArrayView<dimensions, Scalar> m_view;
};


