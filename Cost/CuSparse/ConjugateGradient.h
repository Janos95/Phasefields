//
// Created by janos on 12/7/20.
//

#pragma once

#include <Eigen/SparseCore>

namespace Phasefield {

class ConjugateGradient {

    explicit ConjugateGradient() = default;
    void solve();

private:

    int M = 0, N = 0, nz = 0, *I = NULL, *J = NULL;
    float *val = NULL;
    const float tol = 1e-5f;
    const int max_iter = 10000;
    float *x;
    float *rhs;
    float a, b, na, r0, r1;
    int *d_col, *d_row;
    float *d_val, *d_x, dot;
    float *d_r, *d_p, *d_Ax;
    int k;
    float alpha, beta, alpham1;

};

}


