//
// Created by janos on 08.05.20.
//

#include "normalizeInto.hpp"
#include <limits>
#include <Magnum/Math/FunctionsBatch.h>

using namespace Magnum;
using namespace Corrade;

void normalizeInto(
        Containers::ArrayView<const Double> const& xs,
        Containers::StridedArrayView1D<Float> const& ys,
        Double min, Double max){
    CORRADE_ASSERT(xs.size() == ys.size(), "normalizeInto : arrays dont have the same length",);
    auto l = max - min;
    if(l < std::numeric_limits<Double>::epsilon()){ //handle zero gradient corretly
        for (auto& y : ys) y = .5f;
    } else {
        for (int i = 0; i < xs.size(); ++i) {
            ys[i] = static_cast<Float>((xs[i] - min) / l);
        }
    }
}

void normalizeInto(
        Containers::ArrayView<const Double> const& xs,
        Containers::StridedArrayView1D<Float> const& ys){
    auto [min, max] = Math::minmax(xs);
    normalizeInto(xs, ys, min, max);
}
