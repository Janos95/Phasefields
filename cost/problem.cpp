//
// Created by janos on 28.03.20.
//

#include "problem.hpp"

#include <Eigen/Core>

#include <Corrade/Utility/Assert.h>


std::vector<SumProblem::IndividualCost> SumProblem::computeIndividualCosts(double const* parameters) const {
    std::vector<SumProblem::IndividualCost> costs;
    for(const auto& [name, w, pb, l]: m_problems){
        if(!l)
            continue;
        double cost = 0;
        pb->Evaluate(parameters, &cost, nullptr);
        double out[3];
        l->Evaluate(cost, out);
        auto it = std::lower_bound(costs.begin(), costs.end(),
                [](const auto& f1, const auto& f2){ return f1.name < f2.name; }); //@todo use projection?
        if(it != costs.end() && it->name == name){
            it->cost += cost;
        } else {
            costs.insert(it, IndividualCost{name, w});
        }
    }
    return costs;
}

bool SumProblem::Evaluate(const double* parameters,
                          double* cost,
                          double* jacobians) const {
    auto n = NumParameters();
    *cost = 0.;
    if(jacobians)
        std::fill_n(jacobians, n, 0.);

    auto singleJac = jacobians ? new double[n] : nullptr;

    //@todo this is pretty trivial to thread...
    for(const auto& [_, w, pb, l]: m_problems){
        if(!l)
            continue;
        double residual = 0;
        pb->Evaluate(parameters, &residual, singleJac);
        double out[3];
        l->Evaluate(residual, out);
        if(jacobians){
            Eigen::Map<Eigen::VectorXd> mappedJac(jacobians, n), mappedSingleJac(singleJac, n);
            mappedJac.noalias() += out[1] * mappedSingleJac;
        }
        *cost += out[0];
    }

    delete[] singleJac;

    return true;
}

int SumProblem::NumParameters() const {
    if(m_problems.empty())
        return 0;
    auto numParams = m_problems.front().problem->NumParameters();
#ifndef NDEBUG
    for(const auto& [name, w, prob, l] : m_problems)
        CORRADE_INTERNAL_ASSERT(prob->NumParameters() == numParams);
#endif
    return numParams;
}
