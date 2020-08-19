//
// Created by janos on 08.06.20.
//

#pragma once

#include <type_traits>

template<class F>
struct YCombinator {
    F&& f;

    template<class... Args>
    decltype(auto) operator()(Args&& ... args) const {
        return f(*this, static_cast<Args&&>(args)...);
    }
};

/* deduction guide */
template<class F>
YCombinator(F&& f) -> YCombinator<std::decay_t<F>>;