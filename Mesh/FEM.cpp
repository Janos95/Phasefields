//
// Created by janos on 9/9/20.
//

#include "FEM.h"
#include "SparseMatrix.h"
#include "Mesh.h"

namespace Phasefield {

void IntegralOperatorFeature::update() {
    Mesh& m = mesh();

    arrayResize(m.integral, m.vertexCount());

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

    arrayResize(m.gradient, m.halfEdgeCount());

    for(Face face : m.faces()) {
        Vector3d normal = face.normal();
        double dblA = face.area();
        for(HalfEdge he : face.halfEdges()) {
            Vector3d v = he.asVector();
            m.gradient[he] = Math::cross(normal, v).normalized()*v.length()/dblA;
        }
    }
}

//void StiffnessMatrixFeature::update() {
//
//    auto m = triangles.size();
//    auto n = vertices.size();
//    Array<Triplet> triplets;
//    arrayReserve(triplets, 12*m);
//
//    // run over all faces
//    for(std::size_t i = 0; i < m; ++i) {
//        auto const& t = triangles[i];
//
//        auto a = vertices[t[1]] - vertices[t[2]];
//        auto b = vertices[t[2]] - vertices[t[0]];
//        auto c = vertices[t[0]] - vertices[t[1]];
//
//        auto area = Math::cross(a, -b).length()*.5;
//        Vector3d lengthSqr(a.dot(), b.dot(), c.dot());
//
//        // add local contribution
//        for(int j = 0; j < 3; j++) {
//            auto entry = 0.125*(lengthSqr[j] - lengthSqr[(j + 1)%3] - lengthSqr[(j + 2)%3])/area;
//            arrayAppend(triplets, InPlaceInit, t[(j + 1)%3], t[(j + 2)%3], entry);
//            arrayAppend(triplets, InPlaceInit, t[(j + 2)%3], t[(j + 1)%3], entry);
//            arrayAppend(triplets, InPlaceInit, t[(j + 1)%3], t[(j + 1)%3], -entry);
//            arrayAppend(triplets, InPlaceInit, t[(j + 2)%3], t[(j + 2)%3], -entry);
//        }
//    }
//
//    return triplets;
//}

}