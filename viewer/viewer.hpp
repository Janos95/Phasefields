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

using namespace Math::Literals;

class Viewer: public Magnum::Platform::Application{
public:
    explicit Viewer(int argc, char** argv);

    struct AbstractEventHandler{
        virtual void tickEvent(Scene& scene){ }
        virtual void drawImGui(){ }
        virtual void viewportEvent(ViewportEvent& event, Viewer&){ }
        virtual void keyPressEvent(KeyEvent& event, Viewer&){ }
        virtual void mousePressEvent(MouseEvent& event, Viewer&){ }
        virtual void mouseReleaseEvent(MouseEvent& event, Viewer&){ }
        virtual void mouseMoveEvent(MouseMoveEvent& event, Viewer&){}
        virtual void mouseScrollEvent(MouseScrollEvent& event, Viewer&){ }
        virtual void keyReleaseEvent(KeyEvent& event, Viewer&){ }
        virtual void textInputEvent(TextInputEvent& event, Viewer&){ }
    };

    template<class... H>
    auto insertEventCallbacks(H&&... handlers) {
        (m_eventCallbacks.emplace_back(new(new char[sizeof(H)]) H(std::move(handlers))), ...);
    }

    Viewer& init();
    ArcBallCamera& camera() { return *m_camera; }
    void setScene(Scene& scene) { m_scene = std::addressof(scene); }

private:

    void drawEvent() override;
    void tickEvent() override;
    void viewportEvent(ViewportEvent& event) override;
    void keyPressEvent(KeyEvent& event) override;
    void mousePressEvent(MouseEvent& event) override;
    void mouseReleaseEvent(MouseEvent& event) override;
    void mouseMoveEvent(MouseMoveEvent& event) override;
    void mouseScrollEvent(MouseScrollEvent& event) override;
    void keyReleaseEvent(KeyEvent& event) override;
    void textInputEvent(TextInputEvent& event) override;

    std::vector<std::unique_ptr<AbstractEventHandler>> m_eventCallbacks;

    Corrade::Containers::Optional<ArcBallCamera> m_camera;
    Scene* m_scene = nullptr;
    ImGuiIntegration::Context m_imgui{NoCreate};
    bool m_trackingMouse = false;

};