//
// Created by janos on 11/2/20.
//

#pragma once

#include "Traits.h"

namespace Phasefield {

template<class T1, class T2>
struct Pair {
    T1 first;
    T2 second;
};

template<class T>
void pf_swap(T& a, T& b) {
    T tmp = MOVE(a);
    a = MOVE(b);
    b = MOVE(tmp);
}

}

