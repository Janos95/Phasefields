//
// Created by janos on 05.03.20.
//

#include <Eigen/Core>

double triangleDiameter(
        const Eigen::Vector3d& a,
        const Eigen::Vector3d& b,
        const Eigen::Vector3d& c)
{
    double e1 = (a-b).squaredNorm();
    double e2 = (c-b).squaredNorm();
    double e3 = (a-c).squaredNorm();
    return std::sqrt(std::max({e1, e2, e3}));
}
