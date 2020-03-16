//
// Created by janos on 09.03.20.
//

#include "phasefield_initialization.hpp"

#include <igl/doublearea.h>
#include <random>

auto initAlongZAxis(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::VectorXd& U)
{
    Eigen::VectorXd dblA;
    igl::doublearea(V, F, dblA);
    auto area = dblA.sum();

    // Find the bounding box
    Eigen::Vector3d m = V.colwise().minCoeff();
    Eigen::Vector3d M = V.colwise().maxCoeff();

    double l = m[2], u = M[2], z;

    for (int j = 0; j < 10; ++j) {

        z = .5 * (u + l);

        double s = 0;
        for(int i = 0; i < F.rows(); ++i){
            auto f = F.row(i);
            auto max = std::max({V.row(f[0])[2], V.row(f[1])[2], V.row(f[2])[2]});
            if(max < z)
                s += dblA[i];
        }

        if(s > area / 2.)
            u = z;
        else if(s < area / 2.)
            l = z;
    }

    U = Eigen::VectorXd(V.rows());
    l = m[2];
    u = M[2];

    for(int i = 0; i < V.rows(); ++i){
        double vz = V.row(i)[2];
        U[i] = (vz - z) / (u - l);
    }

}

void initRandomNormal(Eigen::VectorXd& V)
{
    std::default_random_engine engine(0);
    std::normal_distribution distr(.0, .1);

    for(auto& v: V)
        v = distr(engine);
}
