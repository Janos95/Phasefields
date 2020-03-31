//
// Created by janos on 28.03.20.
//

#pragma once


#include <ceres/first_order_function.h>
#include <ceres/loss_function.h>
#include <ceres/iteration_callback.h>

#include <memory>
#include <vector>
#include <folly/Function.h>

struct IterationCallbackWrapper : public ceres::IterationCallback{

    template<class F>
    IterationCallbackWrapper(F&& f) : callback((F&&)f){}

    ceres::CallbackReturnType operator()(ceres::IterationSummary const& summary) override {
         return callback(summary);
     }
     folly::Function<ceres::CallbackReturnType(ceres::IterationSummary const&)> callback;
};

class SumProblem : public ceres::FirstOrderFunction{
public:

    template<class T, class... Args>
    void emplace_back(std::string name, double weight, Args&&... args){
        auto up = std::make_unique<T>((Args&&)args...);
        std::unique_ptr<ceres::ScaledLoss> ul = nullptr;
        if(std::abs(weight) > std::numeric_limits<double>::epsilon()){
            ul = std::make_unique<ceres::ScaledLoss>(new ceres::TrivialLoss, weight, ceres::TAKE_OWNERSHIP);
        }
        m_problems.emplace_back(Functional{std::move(name), weight, std::move(up), std::move(ul)});
    }

    void setWeight(std::string_view name, double weight);

    template<class F>
    void visit(std::string_view name, F f){
        for(auto& [nameOther, _1, p, _2] : m_problems){
            if(nameOther == name){
                f(p);
            }
        }
    }

    void clear() { m_problems.clear(); }

    struct IndividualCost
    {
        std::string_view name;
        double cost;
    };

    [[nodiscard]] std::vector<IndividualCost> computeIndividualCosts(double const* parameters) const;

    bool Evaluate(const double* parameters,
                  double* cost,
                  double* jacobians) const override;

    [[nodiscard]] int NumParameters() const override;

private:
    struct Functional{
        std::string name;
        double weight;
        std::unique_ptr<ceres::FirstOrderFunction> problem;
        std::unique_ptr<ceres::LossFunction> loss = nullptr;
    };

    std::vector<Functional> m_problems;
};