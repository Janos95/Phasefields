//
// Created by janos on 2/23/20.
//

#pragma once

#include "Types.h"
#include "Surface.h"

#include <Corrade/Utility/Assert.h>
#include <Corrade/Containers/Array.h>

namespace Cr = Corrade;

namespace Phasefield {

struct StoppingCriteria {
    StoppingCriteria() = default;

    StoppingCriteria(Face source, size_t numComponents, FaceData<size_t>& components);

    bool operator()(Face target);

    void checkIfFoundAll() const;

    [[nodiscard]] Face target(size_t i) const;

    size_t m_startComponent = 0;
    size_t m_componentsFound = 0;
    FaceData<size_t>* m_components = nullptr;

    Array<bool> m_found;
    Array<Face> m_targetVertices;
};

}


