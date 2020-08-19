//
// Created by janos on 8/17/20.
//

#ifndef PHASEFIELD_FASTMARCHINGMETHOD_H
#define PHASEFIELD_FASTMARCHINGMETHOD_H

#include "Mesh.h"

#include <Corrade/Containers/Containers.h>
#include <Corrade/Containers/Pointer.h>
#include <Magnum/Magnum.h>

namespace Phasefield {

namespace Mg = Magnum;
namespace Cr = Corrade;

class FastMarchingMethod {
public:
    explicit FastMarchingMethod(Mesh&);

    void setSource(Mg::UnsignedInt idx);

    bool step(Mg::UnsignedInt& idx, Mg::Double& distance);

    void reset();

private:

    struct Impl;
    Cr::Containers::Pointer<Impl> m_impl;
};

}

#endif //PHASEFIELD_FASTMARCHINGMETHOD_H
