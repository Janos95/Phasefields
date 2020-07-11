//
// Created by janos on 6/25/20.
//

#pragma once

#include "types.hpp"

#include <Corrade/Containers/Array.h>

struct SparseMatrix {

    std::size_t numRows, numCols;

    Ctrs::Array<double> values;
    Ctrs::Array<std::size_t> rows;
    Ctrs::Array<std::size_t> cols;

    std::size_t size() { return values.size(); }

    bool isDense = false;


    struct RowRange{
        double* b, *e;
        double* begin() { return b; }
        double* end() { return e; }
    };

    RowRange row(std::size_t r){
        auto b = std::lower_bound(rows.begin(), rows.end(), r);
        auto e = std::upper_bound(b, rows.end(), r);
        return {values.begin() + (b - rows.begin()), values.begin() + (e - rows.begin())};
    }


    Ctrs::Array<double> reduceRowwise(){
         Ctrs::Array<double> rowsum(Ctrs::ValueInit, numCols);
         for (int i = 0; i < size(); ++i) {
             rowsum[cols[i]] += values[i];
         }
         return rowsum;
    }
};
