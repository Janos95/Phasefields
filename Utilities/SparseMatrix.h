//
// Created by janos on 6/25/20.
//

#pragma once

#include "Enums.h"
#include "Types.h"

#include <Corrade/Containers/Array.h>
#include <Corrade/Containers/StridedArrayView.h>
#include <Magnum/Types.h>

namespace Phasefield {

namespace Mg = Magnum;

struct Triplet {
    size_t row, column;
    double value;
};

struct SparseMatrix {

    SparseMatrix() = default;

    explicit SparseMatrix(Array<Triplet> triplets);

    size_t numRows, numCols, nnz;

    Array<double> values;
    Array<size_t> rows;
    Array<size_t> cols;

    struct RowRange {
        int current, rowEnd;
    };

    ArrayView<double> row(std::size_t r);

    void clear();

    Array<double> reduceRowwise();
};

}
