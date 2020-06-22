//
// Created by janos on 17.05.20.
//

#pragma once


#include <Magnum/Magnum.h>
#include <Corrade/Containers/StridedArrayView.h>

namespace Mg = Magnum;
namespace Cr = Corrade;

namespace {

    template<class T, class U = void>
    struct choose_tag {
        using type = Cr::Containers::ValueInitT;
    };

    template<class T>
    struct choose_tag<T, std::enable_if_t<std::is_trivial_v<T>>> {
        using type = Cr::Containers::NoInitT;
    };

template<class T>
using choose_tag_t = typename choose_tag<T>::type;

template<unsigned dimensions, class Scalar>
using SizeType = Cr::Containers::StridedDimensions<dimensions, Scalar>;

template<unsigned dim, class Scalar>
constexpr unsigned length(SizeType<dim, Scalar> const& size){
    unsigned l = 1;
    for (auto s : size) l *= s;
    return l;
}

template<
        unsigned... I,
        unsigned dimensions,
        class Scalar,
        class... Idx>
constexpr unsigned index_into(
        std::index_sequence<dimensions>,
        SizeType<dimensions, Scalar> const& strides,
        Idx... idx){
    return  (0 + ... + (idx * strides[I]));
}

}

