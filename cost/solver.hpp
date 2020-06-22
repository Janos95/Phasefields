// int reated by janos on 20.04.20.
//

#pragma once

#include "smart_enum.hpp"
#include "function_ref.hpp"

#include <Corrade/Containers/Array.h>
#include <Magnum/Magnum.h>

namespace Cr = Corrade;
namespace Mg = Magnum;

namespace solver{

    class Problem;

    SMART_ENUM(Solver, Mg::UnsignedInt ,
               IPOPT,
               CERES)

    SMART_ENUM(Status, Mg::UnsignedInt ,
               CONTINUE,
               ABORT,
               FINISHED)

    SMART_ENUM(LineSearchDirection, Mg::UnsignedInt ,
               STEEPEST_DESCENT,
               NONLINEAR_CONJUGATE_GRADIENT,
               LBFGS,
               BFGS,
               SR1)

    struct IterationSummary {
    };

    using iteration_callback_type = function_ref<solver::Status::Value(solver::IterationSummary const&)>;

    struct Options {
        int max_num_iterations = 100;
        bool minimizer_progress_to_stdout = true;
        bool update_state_every_iteration = true;
        LineSearchDirection::Value line_search_direction = LineSearchDirection::LBFGS;
        Solver::Value solver = Solver::CERES;
        Cr::Containers::Array<iteration_callback_type> callbacks;
    };

    struct Summary {

    };

}

void solve(solver::Options& options, solver::Problem& problem, double*, solver::Summary* summary = nullptr);


