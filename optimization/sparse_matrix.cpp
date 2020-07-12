//
// Created by janos on 7/12/20.
//

#include "sparse_matrix.h"

#include <algorithm>

SparseMatrix::RowRange SparseMatrix::row(std::size_t r){
    auto b = std::lower_bound(rows.begin(), rows.end(), r);
    auto e = std::upper_bound(b, rows.end(), r);
    return {values.begin() + (b - rows.begin()), values.begin() + (e - rows.begin())};
}


Containers::Array<double> SparseMatrix::reduceRowwise(){
    Containers::Array<double> rowsum(Containers::ValueInit, numCols);
    for (int i = 0; i < size(); ++i) {
        rowsum[cols[i]] += values[i];
    }
    return rowsum;
}
