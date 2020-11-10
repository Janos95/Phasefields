//
// Created by janos on 9/9/20.
//

#include "FEM.h"
#include "SparseMatrix.h"
#include "Mesh.h"

#include "FEM.hpp"

namespace Phasefield {

void IntegralOperatorFeature::update() {
    Mesh& mesh = getMesh();

    mesh.integral = VertexData<double>{mesh.vertexCount()};

    for(Face face : mesh.faces()) {
        double area = face.area();
        for(Vertex v : face.vertices())
            mesh.integral[v] += area/3.;
    }
}

void MassMatrixFeature::update() {
    //Mesh& mesh = getMesh();

    //Eigen::VectorXd integral{mesh.vertexCount()};
    //integral.setZero();
    //for(Face face : mesh.faces()) {
    //    double area = face.area();
    //    for(Vertex v : face.vertices()) {
    //        integral[v.idx] += area/3.;
    //    }
    //}

    //if(!mesh.fem) {
    //    mesh.fem = pointer<FEM>();
    //}

    //mesh.fem->mass = integral.asDiagonal();
}

void GradientFeature::update() {
    Mesh& mesh = getMesh();

    mesh.gradient = HalfEdgeData<Vector3d>{mesh.halfEdgeCount()};

    for(Face face : mesh.faces()) {
        Vector3d normal = face.normal();
        double dblA = 2*face.area();
        for(HalfEdge he : face.halfEdges()) {
            Vector3d v = he.asVector();
            mesh.gradient[he] = Math::cross(normal, v).normalized()*v.length()/dblA;
        }
    }
}

void StiffnessMatrixFeature::update() {
    Mesh& mesh = getMesh();

    mesh.stiffnessElements = HalfEdgeData<double>{mesh.halfEdgeCount()};
    auto& elements = mesh.stiffnessElements;

    for(HalfEdge he : mesh.halfEdges()) {
        if(he.onBoundaryLoop()) continue;
        HalfEdge he1 = he.next();
        HalfEdge he2 = he1.next();
        elements[he] = Math::dot(mesh.gradient[he1], mesh.gradient[he2])*he.face().area();
    }

    //Array<Eigen::Triplet<double>> triplets{12*mesh.faceCount()};

    //size_t i = 0;
    //for(Face face : mesh.faces()) {
    //    for(HalfEdge he1 : face.halfEdges()) {
    //        Vertex v = he1.next().tip();
    //        for(HalfEdge he2 : face.halfEdges()) {
    //            Vertex w = he2.next().tip();
    //            double value = Math::dot(mesh.gradient[he2], mesh.gradient[he1])*face.area();
    //            triplets[i++] = Eigen::Triplet<double>{v.idx, w.idx, value};
    //        }
    //    }
    //}

    //Eigen::SparseMatrix<double> stiffness(n,n);
    //stiffness.setFromTriplets(triplets.begin(), triplets.end());

    //if(!mesh.fem) {
    //    mesh.fem = pointer<FEM>();
    //}
    //mesh.fem->stiffness = std::move(stiffness);
}

}