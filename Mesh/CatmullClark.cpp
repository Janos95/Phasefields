//
// Created by janos on 9/22/20.
//

#include "CatmullClark.h"

namespace Phasefield {

Mesh catmullClark(Mesh& mesh) {

    FaceData<Vector3> facePoints{mesh.faceCount()};
    EdgeData<Vector3> edgePoints{mesh.edgeCount()};
    VertexData<Vector3> positions{mesh.edgeCount()};

    Array<UnsignedInt> indices;
    auto triangles = arrayCast<Vector3ui>(indices);

    for(Face face : mesh.faces()) {
        Vector3 average{};
        for(Vertex v : face.vertices())
            average += v.position();
        facePoints[face] = average*1./3.;
    }

    for(Edge edge : mesh.edges()) {
        size_t count = 2;
        Vector3 average = edge.vertex2().position() + edge.vertex1().position();
        for(Face f : {edge.halfEdge().face(), edge.halfEdge().twin().face()}) {
            if(f) {
                average += facePoints[f];
                ++count;
            }
        }

        edgePoints[edge] = average*1./float(count);
    }

    for(Vertex vertex : mesh.vertices()) {
        size_t n = float(vertex.degree());
        Vector3 F, R;
        for(HalfEdge he : vertex.outgoingHalfEdges()) {
            Edge e = he.edge();
            R += (e.vertex2().position() + e.vertex1().position())*0.5;
            F += facePoints[he.face()];
        }
        R /= n;
        F /= n;

        positions[vertex] = (F + 2*R + (n - 3)*vertex.position())/n;
    }

    Array<std::pair<size_t, size_t>> edges;


    return Mesh{};
}


}