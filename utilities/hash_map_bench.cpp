//
// Created by janos on 22.05.20.
//

#include "hash_map.hpp"
#include "hash.h"
#include <Magnum/Trade/MeshData.h>
#include <Magnum/Magnum.h>
#include <Magnum/Math/Vector3.h>
#include <Magnum/Primitives/UVSphere.h>
#include <unordered_map>
#include <Magnum/Math/Tags.h>

using namespace Magnum;

struct Edge {
    unsigned int a, b;


    bool operator==(const Edge& other) const {
        return a == other.a && b == other.b;
    }
};

struct HalfEdge {
    static constexpr UnsignedInt Invalid = ~0u;

    UnsignedInt next;
    UnsignedInt face;
    UnsignedInt opposite;
    UnsignedInt vertex;
};

struct Face {
    UnsignedInt halfEdge;
    UnsignedInt degree;
};

struct Mesh {

    using half_edge_descriptor = UnsignedInt;
    using vertex_descriptor = UnsignedInt;
    using face_descriptor = UnsignedInt;

    Trade::MeshData meshData;
    Containers::Array<Face> faces;
    Containers::Array<HalfEdge> halfedges;

    Mesh(Trade::MeshData md) :
            meshData(std::move(md)), halfedges(Containers::NoInit, meshData.indexCount()),
            faces(Containers::NoInit, meshData.indexCount()/3) {
        CORRADE_ASSERT(meshData.primitive() == MeshPrimitive::Triangles,
                       "Can only construct mesh from triangle primitive",);

        auto hash = [](Edge const& e) noexcept { return hash_int(reinterpret_cast<uint64_t const&>(e)); };
        std::unordered_map<Edge, int, decltype(hash)> edges(meshData.indexCount(), hash);
        //HashMap<Edge, int, decltype(hash)> edges(meshData.indexCount(), hash);
        auto triangles = Containers::arrayCast<const Vector3ui>(meshData.indexData());
        for(int i = 0; i < triangles.size(); ++i){
            for(int j = 0; j < 3; ++j){
                auto& he = halfedges[3*i + j];
                he.face = i;
                auto k = (j + 1)%3;
                he.next = 3*i + k;
                he.vertex = triangles[i][k];
                he.opposite = HalfEdge::Invalid;
                auto[it, inserted] = edges.try_emplace(Edge{triangles[i][j], triangles[i][k]}, i);
                if(!inserted){
                    auto at = it->second; /* adjacent triangle */
                    for(int l = 0; l < 3; ++l){
                        if(halfedges[3*at + l].next == triangles[i][j]){
                            he.opposite = 3*at + l;
                            halfedges[3*at + l].opposite = 3*i + j;
                            break;
                        }
                    }
                }
            }
        }
    }

    friend vertex_descriptor add_vertex(Mesh& mesh);

    friend half_edge_descriptor add_edge(Mesh& mesh, vertex_descriptor a, vertex_descriptor b);

    void triangulate() {
        for(auto[hed, degree] : faces){
            if(degree == 3) continue;
            Vector3d mid{Math::ZeroInit};
            auto vd = add_vertex(*this);
            for(UnsignedInt i = 0; i < degree; ++i){
                auto const& he = halfedges[hed];
                add_edge(*this, he.vertex, vd);
                //mid += vertices[he.vertex];
                hed = he.next;
            }
            mid *= 1./degree;
        }
    }
};

Mesh::vertex_descriptor add_vertex(Mesh& mesh) {
    //mesh.meshData.
    //Containers::arrayAppend(mesh.vertices, Containers::InPlaceInit, Math::ZeroInit);
    return mesh.meshData.vertexCount() - 1;
}

Mesh::half_edge_descriptor add_edge(Mesh& mesh, Mesh::vertex_descriptor a, Mesh::vertex_descriptor b) {
    //Containers::arrayAppend(mesh.halfedges, Containers::InPlaceInit, next, face, opposite, b);
    return mesh.halfedges.size() - 1;
}

int main() {

    auto hash = [](int a) noexcept { return hash_int(a); };
    HashMap<int, int, decltype(hash)> map;
    for(int i = 0; i < 1000000; ++i){
        map.tryEmplace(Corrade::Containers::InPlaceInit, rand(), rand());
    }
}