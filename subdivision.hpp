//
// Created by janos on 03.04.20.
//

#pragma once

#include "phasefield_data.hpp"
#include "viewer.hpp"

class Subdivision : public Viewer::AbstractEventHandler {

    PhasefieldData* m_phasefieldData;
    int m_numFaces;
    int m_numVertices;

public:

    void drawEvent() override;
    void subdivide(int);
};



