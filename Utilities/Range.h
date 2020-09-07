//
// Created by janos on 9/6/20.
//

#pragma once

namespace Phasefield {

template<class It>
struct Range {
    It b, e;
    It begin() { return b; }
    It end() { return e; }
};

}

