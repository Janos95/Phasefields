#include "SmoothVertexData.h"
#include "Mesh.h"

#include <Eigen/SparseCore>
#include <Eigen/CholmodSupport>

#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>

namespace Phasefield {

void smoothVertexData(VertexData<double>& vertexData, Mesh& mesh, double c) {
    size_t n = mesh.vertexCount();
    size_t m = mesh.faceCount();

    Eigen::MatrixXd V(n, 3);
    Eigen::MatrixXi F(m, 3);

    for(Vertex v : mesh.vertices()) {
        for(size_t i = 0; i < 3; ++i) {
           V(v.idx, i) = v.position()[i];
        }
    }

    for(Face face : mesh.faces()) {
        size_t i = 0;
        for(Vertex v : face.vertices()) {
            F(face.idx, i++) = v.idx;
        }
    }

    Eigen::SparseMatrix<double> M;
    Eigen::SparseMatrix<double> L;
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
    igl::cotmatrix(V, F, L);

    {
        Eigen::Map<Eigen::VectorXd> b{vertexData.data(), int(n)};
        Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>> solver(M - c*L);
        Eigen::Map<Eigen::VectorXd>{vertexData.data(), int(n)} = solver.solve(b);
    }

    return;

    mesh.requireMassMatrix();
    mesh.requireStiffnessMatrix();

    auto& stiffness = mesh.stiffnessElements;
    auto& mass = mesh.integral;

    Array<Eigen::Triplet<double>> triplets;
    arrayReserve(triplets, mesh.vertexCount() + mesh.halfEdgeCount());

    //for(HalfEdge he : mesh.halfEdges()) {
    //    if(he.edge().onBoundaryLoop()) continue;
    //    arrayAppend(triplets, InPlaceInit, he.tail().idx, he.tip().idx, stiffness[he]);
    //}

    Eigen::VectorXd b(n);

    for(Vertex v : mesh.vertices()) {
        b[v.idx] = vertexData[v];
        arrayAppend(triplets, InPlaceInit, v.idx, v.idx, c*mass[v]);
    }

    Eigen::SparseMatrix<double> A(n,n);
    A.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>> solver(A);
    Eigen::Map<Eigen::VectorXd>{vertexData.data(), int(n)} = solver.solve(b);
}

}