//
// Created by janos on 07.04.20.
//

#pragma once


#include "phasefield_data.hpp"
#include "viewer.hpp"

class ShaderOptions : public Viewer::AbstractEventHandler {
    PhasefieldData& m_phasefieldData;


    Object3D* m_node;



public:

    explicit ShaderOptions(PhasefieldData& data);

    void drawImGui(Viewer&) override;
    void tickEvent(Viewer&) override;
    bool displayPhongOptions();
    bool displayFlatOptions();
    bool displayMeshVisualizerOptions();

};



