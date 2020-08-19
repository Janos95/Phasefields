//
// Created by janos on 22.05.20.
//

#include "HashMap.h"

#include <Corrade/Utility/MurmurHash2.h>

#include <Magnum/Trade/MeshData.h>
#include <Magnum/Magnum.h>
#include <Magnum/Math/Vector3.h>
#include <Magnum/Primitives/UVSphere.h>
#include <Magnum/Math/Tags.h>

using namespace Magnum;

struct Edge {
    UnsignedInt a, b;

    bool operator==(const Edge& other) const {
        return a == other.a && b == other.b;
    }

    std::size_t hash() const {
        const UnsignedLong ab = UnsignedLong(a) + (UnsignedLong(b) << 32ul);
        return *reinterpret_cast<std::size_t const*>(Utility::MurmurHash2{}(reinterpret_cast<const char*>(&ab), sizeof(UnsignedLong)).byteArray());
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

    enum class Attribute : UnsignedShort {
        Position = 1,
        Tangent,
        Bitangent,
        Normal,
        TextureCoordinates,
        Color
    }

    using halfedgedescriptor = unsignedint;
    using vertexdescriptor = unsignedint;
    using facedescriptor = unsignedint;

    Containers::Array<char> vertexData;
    Containers::Array<char> indexData;

    Containers::Array<Trade::MeshAttributeData> attributes;

    Containers::Array<Face> faces;
    Containers::Array<HalfEdge> halfedges;

    Mesh(Trade::MeshData meshData){

        vertexData = meshData.releaseVertexData();
        attributes = meshData.releaseAttributeData();

        auto hash = [](Edge const& e) noexcept { return e.hash(); };
        Containers::HashMap<Edge, int, decltype(hash)> edges(meshData.indexCount());

        //HashMap<Edge, int, decltype(hash)> edges(meshData.indexCount(), hash);
        auto triangles = Containers::arrayCast<const Vector3ui>(indexData);

        for(int i = 0; i < triangles.size(); ++i){
            for(int j = 0; j < 3; ++j){
                auto& he = halfedges[3*i + j];
                he.face = i;
                auto k = (j + 1)%3;
                he.next = 3*i + k;
                he.vertex = triangles[i][k];
                he.opposite = HalfEdge::Invalid;
                auto [it, inserted] = edges.tryEmplace(Edge{triangles[i][j], triangles[i][k]}, i);
                if(!inserted){
                    auto at = it->; /* adjacent triangle */
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

    VertexDescriptor addVertex(){

    }

    HalfEdgeDescriptor addEdge(VertexDescriptor a, VertexDescriptor b){

    }


    void requireEdgeLengths(){

    }

    void requireCornerAngles(){

    }

};

int main() {


}