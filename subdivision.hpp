//
// Created by janos on 03.04.20.
//

#pragma once

#include "phasefield_data.hpp"
#include "viewer.hpp"

class Subdivision : public Viewer::AbstractEventHandler {
    PhasefieldData& m_phasefieldData;
    int m_numFaces = 0;
    int m_numVertices = 0;

    int m_lastNumSubdivision = 0;

public:

    explicit Subdivision(PhasefieldData& data): m_phasefieldData(data) {}
    void drawImGui(Viewer&) override;
    void tickEvent(Viewer&) override {
        if(m_phasefieldData.status == PhasefieldData::Status::NewMesh){
                m_lastNumSubdivision = 0;
                m_numVertices = m_phasefieldData.meshData.vertexCount();
                m_numFaces = m_phasefieldData.meshData.indexCount() / 3;
        }
    }
    void subdivide(int);
};



