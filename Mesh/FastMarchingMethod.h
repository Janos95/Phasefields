//
// Created by janos on 8/17/20.
//

#pragma once

#include "Mesh.h"
#include "MeshElements.h"
#include "Heap.h"

#include <Magnum/Magnum.h>

namespace Phasefield {

namespace Mg = Magnum;
namespace Cr = Corrade;

class FastMarchingMethod {
public:

    explicit FastMarchingMethod(Mesh&);

    void setSource(Vertex v);

    bool step(Vertex& v, Mg::Double& distance);

    void reset();

    void update();

private:

    struct Entry {
        Vertex vertex;
        Double distance;
        auto operator <=>(Entry const& other) const { return distance <=> other.distance; }
    };

    Mesh& m_mesh;
    Heap<Entry> m_frontier;

    VertexData<double> m_distances;
    VertexData<char> m_finalized;

    size_t m_visitedVertexCount;
};

}



