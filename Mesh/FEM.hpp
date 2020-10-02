//
// Created by janos on 10/1/20.
//

#ifndef PHASEFIELD_FEM_HPP
#define PHASEFIELD_FEM_HPP

#include <Eigen/SparseCore>

namespace Phasefield {

struct FEM {
    Eigen::SparseMatrix<double> stiffness;
    Eigen::SparseMatrix<double> mass;
};

}

#endif //PHASEFIELD_FEM_HPP
