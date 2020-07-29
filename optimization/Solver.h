// int reated by janos on 20.04.20.
//

#pragma once

#include "SmartEnum.h"
#include "FunctionRef.hpp"
#include "Optimization.h"

#include <Corrade/Containers/Array.h>
#include <Magnum/Magnum.h>

namespace Cr = Corrade;
namespace Mg = Magnum;

namespace Solver {

SMART_ENUM(Backend, Mg::UnsignedInt,
           IPOPT,
           CERES)

SMART_ENUM(Status, Mg::UnsignedInt,
           CONTINUE,
           ABORT,
           FINISHED)

SMART_ENUM(LineSearchDirection, Mg::UnsignedInt,
           STEEPEST_DESCENT,
           NONLINEAR_CONJUGATE_GRADIENT,
           LBFGS,
           BFGS,
           SR1)

struct IterationSummary {
};

using iteration_callback_type = FunctionRef<Status::Value(IterationSummary const&)>;

struct Options {
    int max_num_iterations = 100;
    bool minimizer_progress_to_stdout = true;
    bool update_state_every_iteration = true;
    LineSearchDirection::Value line_search_direction = LineSearchDirection::LBFGS;
    Backend::Value solver = Backend::CERES;
    Cr::Containers::Array<iteration_callback_type> callbacks;
};

struct Summary {

};

//void solve(Options& options, Problem& problem, double*, Summary* summary = nullptr);

void solve(Options& options, RecursiveProblem& problem, double*, Summary* summary = nullptr);

}



