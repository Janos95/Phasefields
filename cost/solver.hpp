// int reated by janos on 20.04.20.
//

#pragma once

#include "unique_function.h"
#include "problem.hpp"

#include <Corrade/Containers/Array.h>


namespace Cr = Corrade;
namespace Mg = Magnum;

namespace solver{

    enum class Solver : Mg::UnsignedInt {
        IPOPT = 0,
        CERES = 1
    };

    enum class Status : Mg::UnsignedInt {
        CONTINUE = 0,
        ABORT = 1,
        FINISHED = 2
    };

    struct IterationSummary {
    };

    struct Options {
        int max_num_iterations = 100;
        bool minimizer_progress_to_stdout = false;
        bool update_state_every_iteration = true;
        Solver solver = Solver::CERES;
        Cr::Containers::Array<unique_function<Status(IterationSummary const&)>> callbacks;
    };

    struct Summary {

    };

}

void solve(solver::Options& options, solver::Problem& problem, double*, solver::Summary* summary = nullptr);


