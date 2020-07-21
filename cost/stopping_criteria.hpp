//
// Created by janos on 2/23/20.
//

#pragma once


#include <Corrade/Utility/Assert.h>
#include <Corrade/Containers/Array.h>

namespace Cr = Corrade;
namespace Mg = Magnum;

template<class R, class T>
void updateWeight(const int target, const T& w, R& neighbors){
    auto it = neighbors.begin();
    for(; it != neighbors.end(); ++it){
        if(it->vertex == target)
            break;
    }
    CORRADE_INTERNAL_ASSERT(it != neighbors.end());
    it->weight = w;
}

struct StoppingCriteria
{
    StoppingCriteria() = default;
    StoppingCriteria(int source, int numComponents, Cr::Containers::Array<int>& components);

    bool operator()(int target);

    [[nodiscard]] bool foundAll() const;

    [[nodiscard]] int target(int i) const;

    int m_startComponent;
    int m_numComponentsToFind;
    Cr::Containers::Array<int> * m_components;

    int m_numComponentsFound = 0;
    Cr::Containers::Array<bool> m_found;
    Cr::Containers::Array<int> m_targetVertices;
};



