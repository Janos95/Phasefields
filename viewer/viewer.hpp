//
// Created by janos on 08.11.19.
//

#pragma once

#include "arc_ball_camera.hpp"

#include <Corrade/Containers/Reference.h>
#include <Corrade/Containers/Optional.h>
#include <Corrade/Containers/Pointer.h>

#include <Corrade/Containers/GrowableArray.h>

#include <Magnum/Platform/Sdl2Application.h>
#include <Magnum/ImGuiIntegration/Context.h>
#include <Magnum/SceneGraph/RigidMatrixTransformation3D.h>
#include <Magnum/SceneGraph/Object.h>
#include <Magnum/SceneGraph/Scene.h>
#include <Magnum/SceneGraph/Drawable.h>


using namespace Math::Literals;

using Object3D = Magnum::SceneGraph::Object<Magnum::SceneGraph::RigidMatrixTransformation3D>;
using Scene3D = Magnum::SceneGraph::Scene<Magnum::SceneGraph::RigidMatrixTransformation3D>;

class Viewer: public Magnum::Platform::Application{
public:
    explicit Viewer(int argc, char** argv);

    struct AbstractEventHandler{
        virtual void tickEvent(Viewer&){ }
        virtual void drawImGui(Viewer&){ }
        virtual void viewportEvent(ViewportEvent&, Viewer&){ }
        virtual void keyPressEvent(KeyEvent&, Viewer&){ }
        virtual void mousePressEvent(MouseEvent&, Viewer&){ }
        virtual void mouseReleaseEvent(MouseEvent&, Viewer&){ }
        virtual void mouseMoveEvent(MouseMoveEvent&, Viewer&){}
        virtual void mouseScrollEvent(MouseScrollEvent&, Viewer&){ }
        virtual void keyReleaseEvent(KeyEvent&, Viewer&){ }
        virtual void textInputEvent(TextInputEvent&, Viewer&){ }
    };

    auto insertEventCallbacks(std::initializer_list<AbstractEventHandler*> handlers) {
        for(auto handler : handlers)
            Corrade::Containers::arrayAppend(m_eventCallbacks, Containers::InPlaceInit, handler);
    }

    Magnum::SceneGraph::DrawableGroup3D drawableGroup;
    Scene3D scene;
    Corrade::Containers::Optional<ArcBallCamera> camera;

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

    Corrade::Containers::Array<Corrade::Containers::Pointer<AbstractEventHandler>> m_eventCallbacks;

    ImGuiIntegration::Context m_imgui{NoCreate};
    bool m_trackingMouse = false;

};