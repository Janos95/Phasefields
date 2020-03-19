//
// Created by janos on 08.11.19.
//

#pragma once


#include "scene.hpp"
#include "arc_ball_camera.hpp"

#include <Corrade/Containers/Reference.h>
#include <Corrade/Containers/Optional.h>

#include <Magnum/Platform/Sdl2Application.h>
#include <Magnum/ImGuiIntegration/Context.h>

#include <folly/Function.h>

class Viewer: public Magnum::Platform::Application{
public:
    explicit Viewer(int argc, char** argv);

    std::vector<folly::FunctionRef<void(Scene&)>> callbacks;
    Scene scene;

private:

    void drawEvent() override;
    void viewportEvent(ViewportEvent& event) override;
    void keyPressEvent(KeyEvent& event) override;
    void mousePressEvent(MouseEvent& event) override;
    void mouseReleaseEvent(MouseEvent& event) override;
    void mouseMoveEvent(MouseMoveEvent& event) override;
    void mouseScrollEvent(MouseScrollEvent& event) override;

    Corrade::Containers::Optional<ArcBallCamera> m_camera;
    ImGuiIntegration::Context m_imgui{NoCreate};
};