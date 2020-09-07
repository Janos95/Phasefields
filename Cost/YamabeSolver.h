//
// Created by janos on 8/25/20.
//

#ifndef PHASEFIELD_YAMABESOLVER_H
#define PHASEFIELD_YAMABESOLVER_H

#include "../Mesh/Mesh.h"
#include <Ei

namespace Phasefield {

struct YamabeSolver {

    explicit YamabeSolver(Mesh& mesh) : m_mesh(&mesh) {
        m_mesh->requireGaussianCurvature();

    }

    Eigen::SparseMatrix<double> m_A;
    Mesh* m_mesh;
};

}

#endif //PHASEFIELD_YAMABESOLVER_H
