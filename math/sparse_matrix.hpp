//
// Created by janos on 25.05.20.
//

#pragma once

template<class Scalar>
struct Triplet {
    int row, col;
    Scalar scalar;
};

template<class Scalar>
struct SparseMatrix {

    void setFromTriplets(Containers::Array <Triplet<Scalar>> triplets) {
        std::sort(triplets.begin(), triplets.end(),
                  [](auto& x, auto& y) { return std::tie(x.row, x.col) < std::tie(y.row, y.col); });

        Containers::arrayResize(rows, Containers::NoInit, map.size());
        Containers::arrayResize(cols, Containers:NoInit, map.size());
        Containers::arrayResize(values, Containers:NoInit, map.size());
        Containers::arrayResize(values, Containers:NoInit, map.size());
        Containers::arrayResize(starts, Containers::ValueInit, numRows + 1);

        int idx = -1;
        int rowOld = -1, colOld = -1;
        for(auto[row, col, value] : triplets){
            if(row == rowOld && col == colOld){
                values[idx] += value;
                continue;
            }

            ++idx;
            values[idx] = value;
            ++starts[1 + row];

            colOld = column;
            rowOld = row;
        }

        for(int i = 1; i < starts.size(); ++i)
            starts[i] += starts[i - 1];
    }

    int numRows;
    int numCols;

    Containers::Array<int> cols;
    Containers::Array<int> starts;
    Containers::Array <Scalar> values;

};

