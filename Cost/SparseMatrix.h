//
// Created by janos on 6/25/20.
//

#pragma once

#include "Enums.h"

#include <Corrade/Containers/Array.h>
#include <Corrade/Containers/StridedArrayView.h>
#include <Magnum/Types.h>

namespace Phasefield {

namespace Mg = Magnum;

struct Triplet {
    Mg::UnsignedInt row, column;
    double value;
};

struct SparseMatrix {

    SparseMatrix() = default;

    explicit SparseMatrix(Containers::Array<Triplet> triplets);

    int numRows, numCols, nnz;

    Containers::Array<double> values;
    Containers::Array<Mg::UnsignedInt> rows;
    Containers::Array<Mg::UnsignedInt> cols;

    struct RowRange {
        int current, rowEnd;
    };

    Containers::ArrayView<double> row(std::size_t r);

    RowRange rowRange() const;

    void clear();

    Containers::Array<double> reduceRowwise();
};

}
