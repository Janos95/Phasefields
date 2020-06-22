//
// Created by janos on 08.05.20.
//

#pragma once

#include <Corrade/Containers/Containers.h>
#include <Magnum/Magnum.h>

namespace Cr = Corrade;
namespace Mg = Magnum;

void normalizeInto(
        Cr::Containers::ArrayView<const Mg::Double> const& xs,
        Cr::Containers::StridedArrayView1D<Mg::Float> const& ys,
        Mg::Double min, Mg::Double max);


void normalizeInto(
        Cr::Containers::ArrayView<const Mg::Double> const& xs,
        Cr::Containers::StridedArrayView1D<Mg::Float> const& ys);
