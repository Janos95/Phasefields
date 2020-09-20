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

    [[nodiscard]] bool foundAll() const;

    [[nodiscard]] Face target(size_t i) const;

    size_t m_startComponent = 0;
    size_t m_numComponentsToFind = 0;
    FaceData<size_t>* m_components = nullptr;

    size_t m_numComponentsFound = 0;
    Array<char> m_found;
    Array<Face> m_targetVertices;
};

}


