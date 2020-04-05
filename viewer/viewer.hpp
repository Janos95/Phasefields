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
        virtual bool keyPressEvent(KeyEvent& event, Viewer&){ return false; }
        virtual bool mousePressEvent(MouseEvent& event, Viewer&){ return false; }
        virtual bool mouseReleaseEvent(MouseEvent& event, Viewer&){ return false; }
        virtual bool mouseMoveEvent(MouseMoveEvent& event, Viewer&){ return false; }
        virtual bool mouseScrollEvent(MouseScrollEvent& event, Viewer&){ return false; }
        virtual bool keyReleaseEvent(KeyEvent& event, Viewer&){ return false; }
        virtual bool textInputEvent(TextInputEvent& event, Viewer&){ return false; }

        static typename Derived::Status status;
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