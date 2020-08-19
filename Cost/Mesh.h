//
// Created by janos on 8/17/20.
//

#ifndef PHASEFIELD_MESH_H
#define PHASEFIELD_MESH_H

#include <Magnum/Types.h>

namespace Phasefield {

namespace Mg = Magnum;

class Mesh {
public:

    using HalfEdge = Mg::UnsignedInt;
    using Vertex = Mg::UnsignedInt;
    using Face = Mg::UnsignedInt;

    bool isManifold();

    void requireEdgeLengths();

    void requireCornerAngles();

    Mg::UnsignedInt vertexCount();


};

}


#endif //PHASEFIELD_MESH_H
