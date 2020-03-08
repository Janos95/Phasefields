//
// Created by janos on 20.11.19.
//

#include "../cost/quadrature_ref_triangle.hpp"

#include <Eigen/Geometry>
#include <Eigen/Core>

#include <cstdio>

int main()
{

    QuadratureRefTriangle<double> quad;
    auto area1 = quad.integrate([](double w1, double w2){return std::pow(w1 + w2, 3);});
    auto area2 = quad.integrate([&](double w1, double w2){
        auto interp = w1*1 + w2 * 2 + (1. - w1 - w2) * 1;
        return std::pow(interp, 3);});

    printf("quad area is %f\n", area1 + area2);
}