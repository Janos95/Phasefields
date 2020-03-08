//
// Created by janos on 09.12.19.
//

#pragma once

#include <Eigen/Core>
#include <ceres/first_order_function.h>
#include <memory>
#include <Magnum/Magnum.h>
#include <Corrade/Containers/EnumSet.h>


enum class GradientFlag : Magnum::UnsignedByte {
    Analytic = 1 << 0,
    Automatic = 1 << 1,
};

using GradientFlags = Corrade::Containers::EnumSet<GradientFlag>;

CORRADE_ENUMSET_OPERATORS(GradientFlags)

class ConnectednessConstraint : public ceres::FirstOrderFunction
{
public:
    ConnectednessConstraint() noexcept;
    
    ConnectednessConstraint(ConnectednessConstraint&&) noexcept;
    
    ConnectednessConstraint& operator=(ConnectednessConstraint&&);
    
    //deleted copy constructor
    ConnectednessConstraint(ConnectednessConstraint&) = delete;
    
    //deleted copy assigment
    ConnectednessConstraint& operator=(ConnectednessConstraint&) = delete;
    
    ConnectednessConstraint(
            const Eigen::MatrixXd& V,
            const Eigen::MatrixXi& F,
            const double epsilon,
            const double a,
            const double b,
            GradientFlags = GradientFlag::Analytic);

    ~ConnectednessConstraint() override;

    bool Evaluate(double const* parameters, double* cost, double* jacobian) const override ;

    int NumParameters() const override;

private:

    struct ImplAnalytic;
    struct ImplAutodiff;
    std::unique_ptr<ImplAutodiff> m_implAutodiff;
    std::unique_ptr<ImplAnalytic> m_implAnalytic;
};
