// int reated by janos on 20.04.20.
//

#pragma once

#include "function_ref.hpp"

#include <Corrade/Containers/Array.h>
#include <Magnum/Magnum.h>

namespace Cr = Corrade;
namespace Mg = Magnum;

namespace solver{

    class Problem;

    enum class Solver : Mg::UnsignedInt {
        IPOPT = 0,
        CERES = 1
    };

    enum class Status : Mg::UnsignedInt {
        CONTINUE = 0,
        ABORT = 1,
        FINISHED = 2
    };

    enum class LineSearchDirectionType : Mg::UnsignedInt {
        STEEPEST_DESCENT,
        NONLINEAR_CONJUGATE_GRADIENT,
        LBFGS,
        BFGS
    };

    struct IterationSummary {
    };

    using iteration_callback_type = function_ref<solver::Status(solver::IterationSummary const&)>;

    struct Options {
        int max_num_iterations = 100;
        bool minimizer_progress_to_stdout = true;
        bool update_state_every_iteration = true;
        LineSearchDirectionType line_search_direction_type = LineSearchDirectionType::LBFGS;
        Solver solver = Solver::CERES;
        Cr::Containers::Array<iteration_callback_type> callbacks;
    };

    struct Summary {

    };

}

void solve(solver::Options& options, solver::Problem& problem, double*, solver::Summary* summary = nullptr);


