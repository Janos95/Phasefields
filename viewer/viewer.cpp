//
// Created by janos on 08.11.19.
//

#include "viewer.hpp"

#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/Math/Matrix4.h>
#include <Magnum/GL/PixelFormat.h>
#include <Magnum/Math/FunctionsBatch.h>

using namespace Corrade;
using namespace Magnum;

using namespace Math::Literals;


Viewer::Viewer(int argc, char** argv):
    Platform::Application{{argc,argv}, Configuration{}.setWindowFlags(Configuration::WindowFlag::Hidden)
            .setTitle("Viewer")
            .setWindowFlags(Configuration::WindowFlag::Resizable)}
{
    GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
    GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);

    /* Set up the camera */
    {
        /* Setup the arcball after the camera objects */
        const Vector3 eye = Vector3::zAxis(-10.0f);
        const Vector3 center{};
        const Vector3 up = Vector3::yAxis();
        m_camera.emplace(scene.root(), eye, center, up, 45.0_degf,
                               windowSize(), framebufferSize());
    }

    /* Start the timer, loop at 60 Hz max */
    setSwapInterval(1);
    setMinimalLoopPeriod(16);
}


void Viewer::viewportEvent(ViewportEvent& event) {
    GL::defaultFramebuffer.setViewport({{}, event.framebufferSize()});

    m_camera->reshape(event.windowSize(), event.framebufferSize());
    scene.setViewportSize(framebufferSize());
}



void Viewer::keyPressEvent(KeyEvent& event) {
    switch(event.key()) {
        case KeyEvent::Key::L:
            if(m_camera->lagging() > 0.0f) {
                Debug{} << "Lagging disabled";
                m_camera->setLagging(0.0f);
            } else {
                Debug{} << "Lagging enabled";
                m_camera->setLagging(0.85f);
            }
            break;
        case KeyEvent::Key::R:
            m_camera->reset();
            break;

        default: return;
    }

    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
}

void Viewer::mousePressEvent(MouseEvent& event) {
    /* Enable mouse capture so the mouse can drag outside of the window */
    /** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
    SDL_CaptureMouse(SDL_TRUE);

    m_camera->initTransformation(event.position());

    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
}

void Viewer::mouseReleaseEvent(MouseEvent&) {
    /* Disable mouse capture again */
    /** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
    SDL_CaptureMouse(SDL_FALSE);
}

void Viewer::mouseMoveEvent(MouseMoveEvent& event) {
    if(!event.buttons()) return;

    if(event.modifiers() & MouseMoveEvent::Modifier::Shift)
        m_camera->translate(event.position());
    else m_camera->rotate(event.position());

    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
}

void Viewer::mouseScrollEvent(MouseScrollEvent& event) {
    const Float delta = event.offset().y();
    if(Math::abs(delta) < 1.0e-2f) return;

    m_camera->zoom(delta);

    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
}


void Viewer::drawEvent() {
    GL::defaultFramebuffer.clear(GL::FramebufferClear::Color|GL::FramebufferClear::Depth);

    /* Call arcball update in every frame. This will do nothing if the camera
       has not been changed. Otherwise, camera transformation will be
       propagated into the camera objects. */
    bool camChanged = m_camera->update();
    m_camera->draw(scene.drawables());
    swapBuffers();

    if(camChanged) redraw();
}