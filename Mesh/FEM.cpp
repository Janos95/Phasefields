//
// Created by janos on 9/9/20.
//

#include "FEM.h"
#include "SparseMatrix.h"
#include "Mesh.h"

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

    arrayResize(m.massMatrix, 0);
    arrayReserve(m.massMatrix, 9*m.faceCount());

    for(Face face : m.faces()) {
        double area = face.area();

        for(Vertex v : face.vertices()) {
            for(Vertex w : face.vertices()) {
                double factor = v == w  ? 2. : 1.;
                arrayAppend(m.massMatrix, InPlaceInit, v.idx, w.idx, factor*area);
            }
        }
    }
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

    Array<Triplet> triplets;
    arrayReserve(triplets, 12*m.faceCount());

    for(Face face : m.faces()) {

        size_t i = 0;
        size_t t[3];
        Vector3 vertices[3];
        for(Vertex v : face.vertices()) {
            vertices[i] = v.position();
            t[i++] = v.idx;
        }

        Vector3 a = vertices[1] - vertices[2];
        Vector3 b = vertices[2] - vertices[0];
        Vector3 c = vertices[0] - vertices[1];

        double area = double(face.area());
        auto lengthSqr = Vector3d{Vector3(a.dot(), b.dot(), c.dot())};

        // add local contribution
        for(int j = 0; j < 3; j++) {
            double entry = 0.125*(lengthSqr[j] - lengthSqr[(j + 1)%3] - lengthSqr[(j + 2)%3])/area;
            arrayAppend(triplets, InPlaceInit, t[(j + 1)%3], t[(j + 2)%3], entry);
            arrayAppend(triplets, InPlaceInit, t[(j + 2)%3], t[(j + 1)%3], entry);
            arrayAppend(triplets, InPlaceInit, t[(j + 1)%3], t[(j + 1)%3], -entry);
            arrayAppend(triplets, InPlaceInit, t[(j + 2)%3], t[(j + 2)%3], -entry);
        }
    }

    m.stiffnessMatrix = std::move(triplets);
}

}