//
// Created by janos on 9/9/20.
//

#include "MeshFeature.h"
#include "Mesh.h"
#include "MeshElements.h"

#include <Magnum/Math/FunctionsBatch.h>
#include <Corrade/Containers/GrowableArray.h>

namespace Phasefield {

MeshFeature::MeshFeature(Mesh& mesh): m_mesh(mesh) {}

Mesh& MeshFeature::mesh() { return m_mesh; }

void AngleFeature::update() {
    Mesh& m = mesh();
    arrayResize(m.angle, DirectInit, m.halfEdgeCount(), Invalid);
    for(Corner corner : m.corners()) {
        HalfEdge side1 = corner.side1();
        HalfEdge side2 = corner.side2();
        m.angle[corner] = Radd{Mg::Math::angle(side1.direction().normalized(), side2.direction().normalized())};
    }
}

void FaceInformationFeature::update() {
    Mesh& m = mesh();
    arrayResize(m.faceArea, m.faceCount());
    arrayResize(m.faceNormal, m.faceCount());
    arrayResize(m.faceDiameter, m.faceCount());

    m.surfaceArea = 0;
    for(Face face : m.faces()) {
        size_t i = 0;
        Vector3 vertices[3];
        for(Vertex v : face.vertices()) vertices[i++] = v.position();

        Vector3 side1 = vertices[2] - vertices[0];
        Vector3 side2 = vertices[1] - vertices[0];
        Vector3 side3 = vertices[2] - vertices[1];

        Vector3d normal = Vector3d{Math::cross(side1, side2)};

        double area = normal.length()*0.5;
        m.surfaceArea += area;
        m.faceArea[face] = area;
        m.faceNormal[face] = normal.normalized();
        m.faceDiameter[face] = double(Math::max({side1.length(), side2.length(), side3.length()}));
    }
}

void EdgeLengthFeature::update() {
    Mesh& m = mesh();
    arrayResize(m.edgeLength, m.edgeCount());
    for(Edge edge : m.edges())
        m.edgeLength[edge] = double((edge.vertex1().position() - edge.vertex2().position()).length());
}

void GaussianCurvatureFeature::update() {
    Mesh& m = mesh();
    m.requireAngles();
    arrayReserve(m.gaussianCurvature, m.vertexCount());
    for(Vertex vertex : m.vertices()) {
        Radd angleSum{0};
        for(Corner corner : vertex.corners())
            angleSum += corner.angle();
        m.gaussianCurvature[vertex] = 2.*Mg::Math::Constants<double>::pi() - double(angleSum);
    }
}

}


