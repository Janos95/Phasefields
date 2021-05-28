//
// Created by janos on 11/2/20.
//

#pragma once

#include "Traits.h"

namespace Phasefield {

template<class InputIt, class T>
constexpr InputIt find(InputIt first, InputIt last, const T& value) {
    for(; first != last; ++first) {
        if(*first == value) {
            return first;
        }
    }
    return last;
}

template<class InputIt, class UnaryPredicate>
constexpr InputIt findIf(InputIt first, InputIt last, UnaryPredicate p) {
    for(; first != last; ++first) {
        if(p(*first)) {
            return first;
        }
    }
    return last;
}

template<class ForwardIt, class T>
ForwardIt remove(ForwardIt first, ForwardIt last, const T& value) {
    first = find(first, last, value);
    if(first != last)
        for(ForwardIt i = first; ++i != last;)
            if(!(*i == value))
                *first++ = MOVE(*i);
    return first;
}

template<class ForwardIt, class UnaryPredicate>
ForwardIt removeIf(ForwardIt first, ForwardIt last, UnaryPredicate p) {
    first = findIf(first, last, p);
    if(first != last)
        for(ForwardIt i = first; ++i != last;)
            if(!p(*i))
                *first++ = MOVE(*i);
    return first;
}

template<class Rng, class UnaryPredicate>
auto removeIf(Rng&& rng, UnaryPredicate p) {
    return removeIf(rng.begin(), rng.end(), p);
}

template<class Rng, class T>
auto removeIf(Rng&& rng, const T& value) {
    return removeIf(rng.begin(), rng.end(), value);
}

}