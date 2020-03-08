//
// Created by janos on 07.03.20.
//


#pragma once

#include <enoki/autodiff.h>
#include <concepts>

namespace detail {

    template<class T>
    concept IsDiffArray = requires(T x){
        { enoki::detach(x) } -> std::convertible_to<typename T::UnderlyingType>;
        { enoki::gradient(x) } -> std::convertible_to<typename T::UnderlyingType>;
        { enoki::set_requires_gradient(x) };
        { T::simplify_graph_() };
    };

template<IsDiffArray T>
    decltype(auto) detach(T const& arr)
    {
        return enoki::detach(arr);
    }

    //noop
    template<class T>
    auto detach(T x){
        return x;
    }
}

