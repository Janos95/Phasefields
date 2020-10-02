//
// Created by janos on 9/9/20.
//

#include "FEM.h"
#include "SparseMatrix.h"
#include "Mesh.h"

#include "FEM.hpp"

namespace Phasefield {

void IntegralOperatorFeature::update() {
    Mesh& m = mesh();

    m.integral = VertexData<double>{m.vertexCount()};

    for(Face face : m.faces()) {
        double area = face.area();
        for(Vertex v : face.vertices())
            m.integral[v] += area/3.;
    }
}

void MassMatrixFeature::update() {
    Mesh& m = mesh();

    Eigen::VectorXd integral{m.vertexCount()};
    integral.setZero();
    for(Face face : m.faces()) {
        double area = face.area();
        for(Vertex v : face.vertices()) {
            integral[v.idx] += area/3.;
        }
    }

    m.fem->mass = integral.asDiagonal();
}

void GradientFeature::update() {
    Mesh& m = mesh();

    m.gradient = HalfEdgeData<Vector3d>{m.halfEdgeCount()};

    for(Face face : m.faces()) {
        Vector3d normal = face.normal();
        double dblA = 2*face.area();
        for(HalfEdge he : face.halfEdges()) {
            Vector3d v = he.asVector();
            m.gradient[he] = Math::cross(normal, v).normalized()*v.length()/dblA;
        }
    }
}

void StiffnessMatrixFeature::update() {
    Mesh& m = mesh();

    size_t n = m.vertexCount();

    Array<Eigen::Triplet<double>> triplets{12*m.faceCount()};

    size_t i = 0;
    for(Face face : m.faces()) {
        for(HalfEdge he1 : face.halfEdges()) {
            Vertex v = he1.next().tip();
            for(HalfEdge he2 : face.halfEdges()) {
                Vertex w = he2.next().tip();
                double value = Math::dot(m.gradient[he2], m.gradient[he1])*face.area();
                triplets[i++] = Eigen::Triplet<double>{v.idx, w.idx, value};
            }
        }
    }

    Eigen::SparseMatrix<double> stiffness(n,n);
    stiffness.setFromTriplets(triplets.begin(), triplets.end());

    if(!m.fem) {
        m.fem = pointer<FEM>();
    }
    m.fem->stiffness = std::move(stiffness);
}

}