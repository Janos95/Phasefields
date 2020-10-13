//
// Created by janos on 22.04.20.
//

#include "StoppingCriteria.h"
#include "Mesh.h"

#include <Corrade/Utility/FormatStl.h>

using namespace Corrade;

namespace Phasefield {

using Cr::Utility::formatString;

StoppingCriteria::StoppingCriteria(Face source, size_t numComponents, FaceData<size_t>& components) :
        m_startComponent(components[source]),
        m_components(&components),
        m_found(DirectInit, numComponents, false),
        m_componentsFound(components[source] + 1),
        m_targetVertices(DirectInit, numComponents, Invalid)
{
    CORRADE_ASSERT(m_startComponent != Invalid, "Stopping Criteria: source component is not valid",);
}

bool StoppingCriteria::operator()(Face target) {
    size_t comp = (*m_components)[target];

    if(comp == Invalid || comp <= m_startComponent || m_found[comp])
        return false;

    m_found[comp] = true;
    m_targetVertices[comp] = target;
    bool stop = ++m_componentsFound == m_found.size();
    return stop;
}

void StoppingCriteria::checkIfFoundAll() const {
    for(size_t i = m_startComponent + 1; i < m_found.size(); ++i) {
        CORRADE_ASSERT(m_found[i], formatString("Dijkstra did not reach component {} (total {} components) starting from {}", i, m_found.size(), m_startComponent).c_str(),);
    }
}

Face StoppingCriteria::target(size_t i) const {
    CORRADE_INTERNAL_ASSERT(i > m_startComponent);
    Face target = m_targetVertices[i];
    CORRADE_INTERNAL_ASSERT(target);
    return target;
}

}
