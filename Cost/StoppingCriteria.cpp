//
// Created by janos on 22.04.20.
//

#include "StoppingCriteria.h"
#include "Mesh.h"

using namespace Corrade;

namespace Phasefield {

StoppingCriteria::StoppingCriteria(Face source, size_t numComponents, FaceData<size_t>& components) :
        m_startComponent(components[source]),
        m_numComponentsToFind(numComponents - m_startComponent - 1),
        m_components(&components),
        m_found(DirectInit, m_numComponentsToFind + 1, 0),
        m_targetVertices(DirectInit, m_numComponentsToFind + 1, Invalid) {
}

bool StoppingCriteria::operator()(Face target) {
    size_t comp = (*m_components)[target];
    size_t idx = comp - m_startComponent;

    if(comp != Invalid || comp <= m_startComponent || m_found[idx])
        return false;

    m_found[idx] = true;
    m_targetVertices[idx] = target;
    bool stop = ++m_numComponentsFound == m_numComponentsToFind;
    return stop;
}

bool StoppingCriteria::foundAll() const {
    for(char found : m_found.slice(1, m_found.size()))
        if(!found) return false;
    return true;
}

Face StoppingCriteria::target(size_t i) const {
    CORRADE_INTERNAL_ASSERT(i > m_startComponent);
    Face target = m_targetVertices[i - m_startComponent];
    CORRADE_INTERNAL_ASSERT(target);
    return target;
}

}
