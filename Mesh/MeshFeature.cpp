//
// Created by janos on 9/9/20.
//

#include "MeshFeature.h"
#include "Mesh.h"
#include "MeshElements.h"

#include "Bfs.h"

#include <Magnum/Math/FunctionsBatch.h>
#include <Corrade/Containers/GrowableArray.h>

namespace Phasefield {

MeshFeature::MeshFeature(Mesh& mesh, bool ownedByMesh): m_mesh(&mesh), m_ownedByMesh(ownedByMesh) {
    arrayAppend(mesh.m_features, this);
}

Mesh& MeshFeature::mesh() { return *m_mesh; }

MeshFeature::~MeshFeature() {
    if(!m_ownedByMesh) { /* if the mesh is now owning the feature we deregister the feature */
        size_t j = Invalid;
        auto& features = m_mesh->m_features;
        for(size_t i = 0; i < features.size(); ++i) {
            if(m_mesh->m_features[i] == this) {
                j = i;
                break;
            }
        }
        if(j != Invalid) {
            std::swap(features[j], features.back());
            arrayResize(features, features.size() - 1);
        }
    }
}

void AngleFeature::update() {
    Mesh& m = mesh();
    m.angle = CornerData<Radd>{m.halfEdgeCount()};
    for(Corner corner : m.corners()) {
        HalfEdge side1 = corner.side1();
        HalfEdge side2 = corner.side2();
        m.angle[corner] = Radd{Mg::Math::angle(side1.direction().normalized(), side2.direction().normalized())};
    }
}

void FaceInformationFeature::update() {
    Mesh& m = mesh();

    m.faceArea = FaceData<double>{m.faceCount()};
    m.faceNormal = FaceData<Vector3d>{m.faceCount()};
    m.faceDiameter = FaceData<double>{m.faceCount()};
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
    m.edgeLength = EdgeData<double>{m.edgeCount()};
    for(Edge edge : m.edges())
        m.edgeLength[edge] = double((edge.vertex1().position() - edge.vertex2().position()).length());
}

void GaussianCurvatureFeature::update() {
    Mesh& m = mesh();
    m.requireAngles();
    m.gaussianCurvature = VertexData<double>{m.vertexCount()};
    m.faceCurvature = FaceData<double>{m.faceCount()};
    CornerData<double> augmentedAngles{m.halfEdgeCount()};
    VertexData<double> angleSums{m.vertexCount()};

    for(Vertex vertex : m.vertices()) {
        double angleSum{0};
        for(Corner corner : vertex.corners()) {
            angleSum += double(corner.angle());
        }
        angleSums[vertex] = angleSum;

        if(!vertex.onBoundary()) {
            m.gaussianCurvature[vertex] = (2.*Mg::Math::Constants<double>::pi() - double(angleSum))/m.integral[vertex];

            for(Corner corner : vertex.corners()) {
                augmentedAngles[corner] = 2*Mg::Math::Constants<double>::pi()*double(corner.angle())/double(angleSum);
            }
        }
    }

    for(Vertex vertex : m.vertices()) {
        if(vertex.onBoundary()) {
            double angleSumAverage = 0;
            size_t count = 0;
            for(Vertex v : vertex.adjacentVertices()) {
                if(!v.onBoundary()) {
                    angleSumAverage += angleSums[v];
                    ++count;
                }
            }

            angleSumAverage /= double(count);
            m.gaussianCurvature[vertex] = (2.*Mg::Math::Constants<double>::pi() - double(angleSumAverage))/m.integral[vertex];

            for(Corner corner : vertex.corners()) {
                augmentedAngles[corner] = 2*Mg::Math::Constants<double>::pi()*double(corner.angle())/double(angleSumAverage);
            }

        }
    }

    double max= 0;
    for(Face f : m.faces()) {
        double deficit = Math::Constants<double>::pi();
        for(Corner c : f.corners()) {
            deficit -= augmentedAngles[c];
        }
        max = Math::max(deficit, max);
        m.faceCurvature[f] = deficit / f.area();
    }
    Debug{} << "Max deficit" << max;
}

}


