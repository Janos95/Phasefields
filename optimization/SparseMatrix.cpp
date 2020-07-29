//
// Created by janos on 7/12/20.
//

#include "SparseMatrix.h"

#include <Corrade/Containers/GrowableArray.h>

#include <algorithm>

Containers::ArrayView<double> SparseMatrix::row(std::size_t r) {
    auto b = std::lower_bound(rows.begin(), rows.end(), r);
    auto e = std::upper_bound(b, rows.end(), r);
    return {values.begin() + (b - rows.begin()), std::size_t(e - b)};
}


Containers::Array<double> SparseMatrix::reduceRowwise() {
    Containers::Array<double> rowsum(Containers::ValueInit, numCols);
    for(int i = 0; i < nnz; ++i){
        rowsum[cols[i]] += values[i];
    }
    return rowsum;
}

void SparseMatrix::clear() {
    nnz = 0;
    Containers::arrayResize(rows, Containers::NoInit, 0);
    Containers::arrayResize(cols, Containers::NoInit, 0);
}