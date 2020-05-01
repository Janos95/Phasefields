//
// Created by janos on 02.04.20.
//

#pragma once

#include "viewer.hpp"
#include "phasefield_data.hpp"

#include "geodesic_algorithm_exact.h"

namespace Cr = Corrade;
namespace Mg = Magnum;

void paint();
void geodesicSearch();

struct Brush : Viewer::AbstractEventHandler {



    explicit Brush(PhasefieldData& data);

    void mousePressEvent(Viewer::MouseEvent& event, Viewer&) override;
    void mouseMoveEvent(Viewer::MouseMoveEvent& event, Viewer&) override;
    void mouseReleaseEvent(Viewer::MouseEvent& event, Viewer&) override;
    void keyPressEvent(Viewer::KeyEvent&, Viewer&) override;
    void keyReleaseEvent(Viewer::KeyEvent&, Viewer&) override;
    void drawImGui(Viewer&) override;
    void tickEvent(Viewer&) override;

    PhasefieldData& pd;

    double phase = 0;
    double targetDist;
    double recursiveFilterFactor = 0.1;
    double distStep = 0.1;
    Mg::Double maxDist = 20.f;
    bool brushing = false;
    bool stop = true;
    std::vector<std::pair<double, int>> distances;
    geodesic::Mesh mesh;
    Vector3d point;
};



