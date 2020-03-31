//
// Created by janos on 08.11.19.
//

#include "viewer.hpp"

#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/Math/Matrix4.h>
#include <Magnum/GL/PixelFormat.h>
#include <Magnum/Math/FunctionsBatch.h>
#include <Magnum/ImGuiIntegration/Context.hpp>

#include <fmt/core.h>


using namespace Corrade;
using namespace Magnum;

using namespace Math::Literals;


Viewer::Viewer(int argc, char** argv):
    Platform::Application{{argc,argv},NoCreate}
{
    /* Setup window */
    {
        const Vector2 dpiScaling = this->dpiScaling({});
        Configuration conf;
        conf.setTitle("Viewer")
                .setSize(conf.size(), dpiScaling)
                .setWindowFlags(Configuration::WindowFlag::Resizable);
        GLConfiguration glConf;
        glConf.setSampleCount(dpiScaling.max() < 2.0f ? 8 : 2);
        if(!tryCreate(conf, glConf)) {
            create(conf, glConf.setSampleCount(0));
        }
    }

    /* Setup ImGui, load a better font */
    {
        ImGui::CreateContext();
        ImGui::StyleColorsDark();

        ImFontConfig fontConfig;
        fontConfig.FontDataOwnedByAtlas = false;
        const Vector2 size = Vector2{windowSize()}/dpiScaling();
        Utility::Resource rs{"fonts"};
        Containers::ArrayView<const char> font = rs.getRaw("SourceSansPro-Regular.ttf");
        ImGui::GetIO().Fonts->AddFontFromMemoryTTF(
                const_cast<char*>(font.data()), Int(font.size()),
                16.0f*framebufferSize().x()/size.x(), &fontConfig);

        m_imgui = ImGuiIntegration::Context{*ImGui::GetCurrentContext(),
                                                  Vector2{windowSize()}/dpiScaling(), windowSize(), framebufferSize()};

        /* Setup proper blending to be used by ImGui */
        GL::Renderer::setBlendEquation(GL::Renderer::BlendEquation::Add,
                                       GL::Renderer::BlendEquation::Add);
        GL::Renderer::setBlendFunction(GL::Renderer::BlendFunction::SourceAlpha,
                                       GL::Renderer::BlendFunction::OneMinusSourceAlpha);

    }

    GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
    GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);

    /* Start the timer, loop at 60 Hz max */
    setSwapInterval(1);
    setMinimalLoopPeriod(16);
}

Viewer& Viewer::init(){
    CORRADE_INTERNAL_ASSERT(m_scene);
    m_scene->setViewportSize(framebufferSize());
    /* Setup the arcball after the camera objects */
    const Vector3 eye = Vector3::zAxis(-10.0f);
    const Vector3 center{};
    const Vector3 up = Vector3::yAxis();
    m_camera.emplace(m_scene->root(), eye, center, up, 45.0_degf,
                     windowSize(), framebufferSize());
    return *this;
}


void Viewer::viewportEvent(ViewportEvent& event) {
    GL::defaultFramebuffer.setViewport({{}, event.framebufferSize()});
    if(m_scene){
        m_camera->reshape(event.windowSize(), event.framebufferSize());
        m_scene->setViewportSize(event.framebufferSize());
    }
    m_imgui.relayout(Vector2{event.windowSize()}/event.dpiScaling(),
                    event.windowSize(), event.framebufferSize());
}


void Viewer::keyPressEvent(KeyEvent& event) {
    if(m_imgui.handleKeyPressEvent(event)) {
        event.setAccepted();
        return;
    }

    if(!m_camera) return;

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

        default:
            if(!m_imgui.handleKeyPressEvent(event))
                return;
            break;
    }

    event.setAccepted();
    redraw(); /* camera or gui element has changed, redraw! */
}

void Viewer::keyReleaseEvent(KeyEvent& event) {
    if(m_imgui.handleKeyReleaseEvent(event)) {
        event.setAccepted();
        return;
    }
}

void Viewer::textInputEvent(TextInputEvent& event) {
    if(m_imgui.handleTextInputEvent(event)) {
        event.setAccepted();
        return;
    }
}

void Viewer::mousePressEvent(MouseEvent& event) {
    if(m_imgui.handleMousePressEvent(event)){
        event.setAccepted();
        return;
    }

    if(!m_camera) return;

    m_trackingMouse = true;
    ///* Enable mouse capture so the mouse can drag outside of the window */
    ///** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
    SDL_CaptureMouse(SDL_TRUE);

    m_camera->initTransformation(event.position());

    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
}

void Viewer::mouseReleaseEvent(MouseEvent& event) {
    if(m_imgui.handleMouseReleaseEvent(event)) {
        event.setAccepted();
        return;
    }

    if(!m_camera) return;
    /* Disable mouse capture again */
    /** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
    if(m_trackingMouse){
        SDL_CaptureMouse(SDL_FALSE);
        m_trackingMouse = false;
        event.setAccepted();
    }
}

void Viewer::mouseMoveEvent(MouseMoveEvent& event) {
    if(m_imgui.handleMouseMoveEvent(event)){
        event.setAccepted();
        return;
    }

    if(!m_camera || !event.buttons()) return;

    if(event.modifiers() & MouseMoveEvent::Modifier::Shift)
        m_camera->translate(event.position());
    else m_camera->rotate(event.position());

    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
}

void Viewer::mouseScrollEvent(MouseScrollEvent& event) {
    if(m_imgui.handleMouseScrollEvent(event)) {
        /* Prevent scrolling the page */
        event.setAccepted();
        return;
    }

    if(!m_camera) return;

    const Float delta = event.offset().y();
    if(Math::abs(delta) < 1.0e-2f) return;

    m_camera->zoom(delta);

    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
}

void Viewer::tickEvent()
{
    for(auto& cb : tickCallbacks)
        cb(m_scene);
    if(m_scene && !m_camera) init();
    if(m_scene && m_scene->isDirty()){
        m_scene->setClean();
        redraw();
    }
}

void Viewer::drawEvent() {
    GL::defaultFramebuffer.clear(GL::FramebufferClear::Color|GL::FramebufferClear::Depth);
    m_imgui.newFrame();

    /* Enable text input, if needed */
    if(ImGui::GetIO().WantTextInput && !isTextInputActive())
        startTextInput();
    else if(!ImGui::GetIO().WantTextInput && isTextInputActive())
        stopTextInput();

    /* draw scene */
    if(m_scene){
        bool camChanged = m_camera->update();
        m_camera->draw(m_scene->drawables());
    }

    for(auto& menuCB : menuCallbacks)
        menuCB(m_imgui);
    m_imgui.updateApplicationCursor(*this);

    /* Render ImGui window */
    {
        GL::Renderer::enable(GL::Renderer::Feature::Blending);
        GL::Renderer::disable(GL::Renderer::Feature::FaceCulling);
        GL::Renderer::disable(GL::Renderer::Feature::DepthTest);
        GL::Renderer::enable(GL::Renderer::Feature::ScissorTest);

        m_imgui.drawFrame();

        GL::Renderer::disable(GL::Renderer::Feature::ScissorTest);
        GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
        GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);
        GL::Renderer::disable(GL::Renderer::Feature::Blending);
    }

    swapBuffers();
    redraw();
}

