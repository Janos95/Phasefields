//
// Created by janos on 07.03.20.
//


#pragma once

namespace enoki {

    template<typename T>
    inline void set_requires_gradient(T &a, bool value);

    template<typename T>
    inline auto detach(T const& a) -> typename T::UnderlyingType;

    template<typename T>
    void backward(const T &a, bool free_graph);

    template<typename T>
    decltype(auto) gradient(T &&value);

    template<class T>
    class DiffArray;
}

#ifdef __cpp_concepts
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

#else
namespace detail
{
    template< class >
    struct is_diff_array : std::false_type { };

    template< class T >
    struct is_diff_array<enoki::DiffArray<T>> : std::true_type { };

    template< class T >
    constexpr bool IsDiffArray = is_diff_array<T>::value;

    template<class T>
    auto detach(T x){
        if constexpr(IsDiffArray<T>)
            return enoki::detach(x);
        else
            return x;
    }
}
#endif


