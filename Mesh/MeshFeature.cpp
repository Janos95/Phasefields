//
// Created by janos on 9/9/20.
//

#include "MeshFeature.h"
#include "Mesh.h"
#include "MeshElements.h"

#include <Corrade/Containers/GrowableArray.h>

namespace Phasefield {

MeshFeature::MeshFeature(Mesh& mesh): m_mesh(mesh) {}

Mesh& MeshFeature::mesh() { return m_mesh; }

void AngleFeature::update() {
    Mesh& m = mesh();
    arrayResize(m.angle, DirectInit, m.angleCount(), Invalid);
    for(Corner corner : m.corners()) {
        HalfEdge side1 = corner.side1();
        HalfEdge side2 = corner.side2();
        m.angle[corner] = Rad{Mg::Math::angle(side1.direction().normalized(), side2.direction().normalized())};
    }
}

void FaceAreaFeature::update() {
    Mesh& m = mesh();
    arrayResize(m.faceArea, m.faceCount());
    arrayResize(m.faceNormal, m.faceCount());

    m.surfaceArea = 0;
    for(Face face : m.faces()) {
        Corner corner = face.halfEdge().corner();
        Vector3 side1 = corner.side1().direction();
        Vector3 side2 = corner.side2().direction();
        Vector3d normal = Vector3d{Math::cross(side1, side2)};
        double area = normal.length()*0.5;
        m.surfaceArea += area;
        m.faceArea[face] = area;
        m.faceNormal[face] = normal.normalized();
    }
}

void EdgeLengthFeature::update() {
    Mesh& m = mesh();
    arrayResize(m.edgeLength, m.edgeCount());
    for(Edge edge : m.edges())
        m.edgeLength[edge] = (edge.vertex1().position() - edge.vertex2().position()).length();
}

void GaussianCurvatureFeature::update() {
    Mesh& m = mesh();
    m.requireAngles();
    arrayReserve(m.gaussianCurvature, m.vertexCount());
    for(Vertex vertex : m.vertices()) {
        Rad angleSum{0};
        for(Corner corner : vertex.corners())
            angleSum += corner.angle();
        m.gaussianCurvature[vertex] = 2.f*Mg::Math::Constants<float>::pi() - static_cast<float>(angleSum);
    }
}

}


