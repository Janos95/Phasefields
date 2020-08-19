//
// Created by janos on 22.04.20.
//

#include "StoppingCriteria.h"

using namespace Corrade;

StoppingCriteria::StoppingCriteria(int source, int numComponents, Cr::Containers::Array<int>& components) :
        m_startComponent(components[source]),
        m_numComponentsToFind(numComponents - m_startComponent - 1),
        m_components(&components),
        m_found(Cr::Containers::DirectInit, m_numComponentsToFind + 1, false),
        m_targetVertices(Cr::Containers::DirectInit, m_numComponentsToFind + 1, -1) {
}

bool StoppingCriteria::operator()(int target) {
    auto comp = (*m_components)[target];
    auto idx = comp - m_startComponent;
    //if target is not in the interface then comp is -1.
    if(comp <= m_startComponent || m_found[idx])
        return false;

    m_found[idx] = true;
    m_targetVertices[idx] = target;
    auto stop = ++m_numComponentsFound == m_numComponentsToFind;
    return stop;
}

bool StoppingCriteria::foundAll() const {
    for(auto it = m_found.begin() + 1; it != m_found.end(); ++it)
        if(!(*it)) return false;
    return true;
}

int StoppingCriteria::target(int i) const {
    CORRADE_INTERNAL_ASSERT(i > m_startComponent);
    auto target = m_targetVertices[i - m_startComponent];
    CORRADE_INTERNAL_ASSERT(target >= 0);
    return target;
}
