//
// Created by janos on 03.04.20.
//

#pragma once

#include "scene.hpp"
#include "viewer.hpp"
#include "phasefield_data.hpp"

struct UpdateScene : public Viewer::AbstractEventHandler {

    PhasefieldData& phasefieldData;

    void tickEvent(Scene& scene) override;

    void makeScene(Scene& scene);

    void reuploadVertices(Scene& scene);

    void uploadIndices(Scene& scene);

};



