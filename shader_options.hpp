//
// Created by janos on 07.04.20.
//

#pragma once


#include "phasefield_data.hpp"
#include "viewer.hpp"

class ShaderOptions : public Viewer::AbstractEventHandler {
    PhasefieldData& m_phasefieldData;
    bool m_updateNode = false;
    DrawableType m_type;

    Float m_lineWidth = 1.f;
    Float m_wireframeWidth = 1.f;
public:

    explicit ShaderOptions(PhasefieldData& data): m_phasefieldData(data) {}

    void drawImGui() override;
    void tickEvent(Scene &scene) override;
    void handlePhong();
    void handleFlat();
    void handleMeshVis();

};



