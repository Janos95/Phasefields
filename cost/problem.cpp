//
// Created by janos on 28.03.20.
//

#include "problem.hpp"

#include <tbb/task_group.h>

#include <Corrade/Utility/Assert.h>
#include <Corrade/Utility/Algorithms.h>

#include <mutex>

using namespace Corrade;
using namespace Magnum;

namespace solver {

bool Problem::evaluate(double const *parameters,
                       double *cost,
                       double *jacobians) const {
    CORRADE_ASSERT(!functionals.empty(), "Problem : no functionals in problem", false);
    const std::size_t n = numParameters();

    *cost = 0.;
    if (jacobians)
        std::fill_n(jacobians, n, 0.);

    std::mutex m1;
    std::mutex m2;

    //tbb::task_group g;
    for (auto const& functional: functionals) {
        //g.run([&] {
            Containers::Array<Double> grad(Containers::NoInit, jacobians ? n : 0);
            Double residual = 0;
            functional->evaluate(parameters, &residual, jacobians ? grad.data() : nullptr);


            Double rho[3], scaling[3] = {residual, 1., 0.};
            if(functional->scaling){
                auto s = *functional->scaling;
                for (auto& r : scaling) r *= s;
            }
            functional->loss->Evaluate(scaling[0], rho);
            if (jacobians) {
                auto derivative = rho[1] * scaling[1];
                std::lock_guard l(m1);
                for (int i = 0; i < n; ++i)
                    jacobians[i] += derivative * grad[i];
            }
            std::lock_guard l(m2);
            *cost += rho[0];
       // });
    }
    //g.wait();
    return true;
}

int Problem::numParameters() const {
    if (functionals.empty() && constraints.empty())
        return 0;
    auto numParams = functionals.front()->numParameters();
#ifndef NDEBUG
    for (auto const& f : functionals)
        CORRADE_INTERNAL_ASSERT(f->numParameters() == numParams);
    for (auto const& c : constraints)
        CORRADE_INTERNAL_ASSERT(c->numParameters() == numParams);
#endif
    return numParams;
}

}
