//
// Created by janos on 11/1/20.
//

#pragma once

#include "Types.h"
#include "RecursiveProblem.h"
#include <Corrade/Containers/ArrayView.h>

namespace Phasefield {

/**
 * This is dumb wrapper around the lbfgs solver (see contrib module).
 * The original solver only allows to call into the solver once
 */
class LbfgsSolver {
public:

    LbfgsSolver(Solver::Options& options, Solver::RecursiveProblem& problem, ArrayView<double> data);

    LbfgsSolver(LbfgsSolver const&) = delete;
    LbfgsSolver(LbfgsSolver&&) = delete;
    LbfgsSolver& operator=(LbfgsSolver&&) = delete;
    LbfgsSolver& operator=(LbfgsSolver const&) = delete;

    ~LbfgsSolver();

    int runOneIteration();

private:

    ArrayView<double> m_x;
    double* m_fx;
    Solver::RecursiveProblem& m_problem;
    Solver::Options& m_options;

    struct Data;
    Data* m_data;

    void* m_cd;

    int ret;
    int i, j, k, ls, end, bound;

    double step;
    double ys, yy;
    double xnorm, gnorm, beta;
    double fx = 0.;
    double rate = 0.;
};

}

