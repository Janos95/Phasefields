//
// Created by janos on 9/9/20.
//

#include "MeshFeature.h"
#include "Mesh.h"
#include "MeshElements.h"

#include <Magnum/Math/FunctionsBatch.h>
#include <Corrade/Containers/GrowableArray.h>

namespace Phasefield {

MeshFeature::MeshFeature(Mesh& mesh, bool ownedByMesh): m_mesh(&mesh), m_ownedByMesh(ownedByMesh) {
    arrayAppend(mesh.m_features, this);
}

Mesh& MeshFeature::getMesh() { return *m_mesh; }

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
    Mesh& mesh = getMesh();
    mesh.angle = CornerData<Radd>{mesh.halfEdgeCount()};
    for(Corner corner : mesh.corners()) {
        HalfEdge side1 = corner.side1();
        HalfEdge side2 = corner.side2();
        mesh.angle[corner] = Radd{Mg::Math::angle(side1.direction().normalized(), side2.direction().normalized())};
    }
}

void FaceInformationFeature::update() {
    Mesh& mesh = getMesh();

    mesh.faceArea = FaceData<double>{mesh.faceCount()};
    mesh.faceNormal = FaceData<Vector3d>{mesh.faceCount()};
    mesh.faceDiameter = FaceData<double>{mesh.faceCount()};
    mesh.surfaceArea = 0;

    for(Face face : mesh.faces()) {
        size_t i = 0;
        Vector3 vertices[3];
        for(Vertex v : face.vertices()) vertices[i++] = v.position();

        Vector3 side1 = vertices[2] - vertices[0];
        Vector3 side2 = vertices[1] - vertices[0];
        Vector3 side3 = vertices[2] - vertices[1];

        Vector3d normal = Vector3d{Math::cross(side1, side2)};

        double area = normal.length()*0.5;
        mesh.surfaceArea += area;
        mesh.faceArea[face] = area;
        mesh.faceNormal[face] = normal.normalized();
        mesh.faceDiameter[face] = double(Math::max({side1.length(), side2.length(), side3.length()}));
    }
}

void EdgeLengthFeature::update() {
    Mesh& mesh = getMesh();
    mesh.edgeLength = EdgeData<double>{mesh.edgeCount()};
    for(Edge edge : mesh.edges())
        mesh.edgeLength[edge] = double((edge.vertex1().position() - edge.vertex2().position()).length());
}

void GaussianCurvatureFeature::update() {
    Mesh& mesh= getMesh();
    mesh.requireAngles();
    mesh.gaussianCurvature = VertexData<double>{mesh.vertexCount()};
    mesh.faceCurvature = FaceData<double>{mesh.faceCount()};
    CornerData<double> augmentedAngles{mesh.halfEdgeCount()};
    VertexData<double> angleSums{mesh.vertexCount()};

    for(Vertex vertex : mesh.vertices()) {
        double angleSum{0};
        for(Corner corner : vertex.corners()) {
            angleSum += double(corner.angle());
        }
        angleSums[vertex] = angleSum;

        if(!vertex.onBoundary()) {
            mesh.gaussianCurvature[vertex] = (2.*Mg::Math::Constants<double>::pi() - double(angleSum))/mesh.integral[vertex];

            for(Corner corner : vertex.corners()) {
                augmentedAngles[corner] = 2*Mg::Math::Constants<double>::pi()*double(corner.angle())/double(angleSum);
            }
        }
    }

    for(Vertex vertex : mesh.vertices()) {
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
            mesh.gaussianCurvature[vertex] = (2.*Mg::Math::Constants<double>::pi() - double(angleSumAverage))/mesh.integral[vertex];

            for(Corner corner : vertex.corners()) {
                augmentedAngles[corner] = 2*Mg::Math::Constants<double>::pi()*double(corner.angle())/double(angleSumAverage);
            }

        }
    }

    double max= 0;
    for(Face f : mesh.faces()) {
        double deficit = Math::Constants<double>::pi();
        for(Corner c : f.corners()) {
            deficit -= augmentedAngles[c];
        }
        max = Math::max(deficit, max);
        mesh.faceCurvature[f] = deficit / f.area();
    }
    Debug{} << "Max deficit" << max;

    //Debug{} << "SETTING GAUSIAN CURVATURE --- ONLY FOR BLENDER TORUS SUITABLE";
    //for(Vertex v : mesh.vertices()) {
    //    auto& p = v.position();
    //    double x = p.x();
    //    double y = p.y();
    //    double z = p.z();

    //    double xy2 = 1.-Math::sqrt(x*x + y*y);
    //    double cosTheta = xy2/Math::sqrt(z*z + xy2*xy2);
    //    mesh.gaussianCurvature[v] = 2.*cosTheta/(1. + cosTheta*0.25);
    //}
}

void BoundaryInformation::update() {
    Mesh& mesh = getMesh();

    arrayResize(mesh.isOnBoundary, NoInit, mesh.vertexCount());

    auto onBoundary = [](Vertex v) {
        for(HalfEdge he : v.outgoingHalfEdges())
            if(he.onBoundaryLoop()) return true;
        return false;
    };

    for(Vertex v : mesh.vertices()) {
        bool isOnBoundary = onBoundary(v);
        mesh.isOnBoundary[v] = isOnBoundary;
    }
}

}


