//
// Created by janos on 6/25/20.
//

#pragma once

#include "types.hpp"

#include <Corrade/Containers/Array.h>
#include <Corrade/Containers/StridedArrayView.h>

struct Triplet {
    std::size_t row, column;
    double value;
};

struct SparseMatrix {

    explicit SparseMatrix(Containers::Array<Triplet> triplets);

    std::size_t numRows, numCols;

    Containers::Array<double> values;
    Containers::Array<std::size_t> rows;
    Containers::Array<std::size_t> cols;

    std::size_t size() { return values.size(); }

    bool isDense = false;


    struct Row{
        double* b, *e;
        double* begin() { return b; }
        double* end() { return e; }

    };

    struct RowRange{
        int current, rowEnd;
    };

    Row row(std::size_t r);
    RowRange rowRange();

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

Vector operator*(SparseMatrix const& m, Vector v){
    Vector result(m.numRows);
    for(auto row : m.rowRange()){
        auto rowIdx = row.idx();
        for(auto colIdx : row){
            result[rowIdx] = m.values[colIdx] * v[m.cols[colIdx]];
        }
    }
    return result;
}
