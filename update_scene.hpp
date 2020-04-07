//
// Created by janos on 03.04.20.
//

#pragma once

#include "scene.hpp"
#include "viewer.hpp"
#include "phasefield_data.hpp"

struct UpdateScene : public Viewer::AbstractEventHandler {

    explicit UpdateScene(PhasefieldData& data): phasefieldData(data) {}

    PhasefieldData& phasefieldData;

    void tickEvent(Scene&scene) override;

    void makeScene(Scene& scene);

    void reuploadVertices();

    void reuploadIndices();

    void reupload();
};



