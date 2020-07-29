//
// Created by janos on 6/25/20.
//

#pragma once

#include "types.hpp"

#include <Corrade/Containers/Array.h>
#include <Corrade/Containers/StridedArrayView.h>

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

struct Matrix {
    Containers::Array<double> values;
    Containers::StridedArrayView2D<double> view;
};

struct Vector {
    explicit Vector(std::size_t size);

    Containers::Array<double> values;
};

struct VectorView {
    Containers::ArrayView<double> view;
};

struct MatrixView {
    Containers::StridedArrayView2D<double> view;
};


