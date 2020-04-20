//
// Created by janos on 28.03.20.
//

#include "problem.hpp"

#include <Eigen/Core>

#include <tbb/task_group.h>

#include <Corrade/Utility/Assert.h>
#include <fmt/core.h>

#include <mutex>

std::vector<SumProblem::IndividualCost> SumProblem::computeIndividualCosts(double const* parameters) const {
    std::vector<SumProblem::IndividualCost> costs;
    for(auto const& [name, w, pb]: functionals){
        if(!l)
            continue;
        double cost = 0;
        pb->Evaluate(parameters, &cost, nullptr);
        double out[3];
        l->Evaluate(cost, out);
        auto it = std::lower_bound(costs.begin(), costs.end(), IndividualCost{name},
                [](auto const& f1, auto const& f2){ return f1.name < f2.name; }); //@todo use projection?
        if(it != costs.end() && it->name == name){
            it->cost += cost;
        } else {
            costs.insert(it, IndividualCost{name, w});
        }
    }
    return costs;
}

void SumProblem::setWeight(std::string_view name, double weight) {
    for(auto& [nameOther, _1, _2, l] : m_problems){
        if(nameOther == name){
           if(std::abs(weight) < std::numeric_limits<double>::epsilon())
               l = nullptr;
           else
               l = std::make_unique<ceres::ScaledLoss>(new ceres::TrivialLoss, weight, ceres::TAKE_OWNERSHIP);
        }
    }
}


bool SumProblem::Evaluate(double const* parameters,
                          double* cost,
                          double* jacobians) const {
    auto n = NumParameters();
    *cost = 0.;
    if(jacobians)
        std::fill_n(jacobians, n, 0.);

    std::mutex mutex;

    tbb::task_group g;
    for(auto const& problem: m_problems){
        if(!problem.loss)
            continue;
        g.run([&] {
            auto singleJac = jacobians ? new double[n] : nullptr;
            double residual = 0;
            problem.problem->Evaluate(parameters, &residual, singleJac);
            double out[3];
            problem.loss->Evaluate(residual, out);
            if (jacobians) {
                Eigen::Map<Eigen::VectorXd> mappedJac(jacobians, n), mappedSingleJac(singleJac, n);
                std::lock_guard l(mutex);
                mappedJac.noalias() += out[1] * mappedSingleJac;
            }
            {
                std::lock_guard l(mutex);
                *cost += out[0];
            }
            delete[] singleJac;
        });
    }
    g.wait();
    return true;
}

int SumProblem::NumParameters() const {
    if(m_problems.empty())
        return 0;
    auto numParams = m_problems.front().problem->NumParameters();
#ifndef NDEBUG
    for(auto const& [name, w, prob, l] : m_problems)
        CORRADE_INTERNAL_ASSERT(prob->NumParameters() == numParams);
#endif
    return numParams;
}
