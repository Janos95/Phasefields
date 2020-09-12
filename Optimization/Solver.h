// int reated by janos on 20.04.20.
//

#pragma once

#include "SmartEnum.h"
#include "UniqueFunction.h"
#include "Optimization.h"
#include "Types.h"

#include <Corrade/Containers/Array.h>
#include <Magnum/Magnum.h>

namespace Cr = Corrade;
namespace Mg = Magnum;

namespace Phasefield::Solver {

SMART_ENUM(Backend, size_t,
           IPOPT,
           CERES)

SMART_ENUM(Status, size_t,
           CONTINUE,
           USER_ABORTED,
           FINISHED)

SMART_ENUM(LineSearchDirection, size_t,
           STEEPEST_DESCENT,
           NONLINEAR_CONJUGATE_GRADIENT,
           LBFGS,
           BFGS,
           SR1)

struct IterationSummary {
};

using IterationCallbackType = UniqueFunction<Status::Value(IterationSummary const&)>;

struct Options {
    size_t max_num_iterations = 100;
    bool minimizer_progress_to_stdout = true;
    bool update_state_every_iteration = true;
    LineSearchDirection::Value line_search_direction = LineSearchDirection::LBFGS;
    Backend::Value solver = Backend::CERES;
    Array<IterationCallbackType> callbacks;
};

struct Summary {

};

//void solve(Options& options, Problem& problem, double*, Summary* summary = nullptr);

void solve(Options& options, RecursiveProblem& problem, ArrayView<double> data, Summary* summary = nullptr);

}



