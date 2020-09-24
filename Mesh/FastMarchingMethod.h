//
// Created by janos on 8/17/20.
//

#pragma once

#include "Mesh.h"
#include "MeshElements.h"
#include "Heap.h"
#include "GraphCommon.h"

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

    Mesh& m_mesh;
    Heap<Graph::HeapElement<Vertex>> m_frontier;

    VertexData<double> m_distances;
    VertexData<char> m_finalized;

};

}



