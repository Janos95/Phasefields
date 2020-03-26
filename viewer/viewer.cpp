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

Viewer& Viewer::setScene(Scene& scene){
    m_scene = std::addressof(scene);

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
    m_camera->reshape(event.windowSize(), event.framebufferSize());
    m_scene->setViewportSize(framebufferSize());
    m_imgui.relayout(Vector2{event.windowSize()}/event.dpiScaling(),
                    event.windowSize(), event.framebufferSize());
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

    /* Disable mouse capture again */
    /** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
    if(m_trackingMouse){
        fmt::print("stopping mouse tracking\n");
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

    if(!event.buttons()) return;

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

    const Float delta = event.offset().y();
    if(Math::abs(delta) < 1.0e-2f) return;

    m_camera->zoom(delta);

    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
}

void Viewer::tickEvent()
{
    for(auto& cb : callbacks)
        cb(*this);
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
    bool camChanged = m_camera->update();
    m_camera->draw(m_scene->drawables());

    //showMenu();
    ImGui::ShowDemoWindow();
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

void Viewer::showMenu() {
    ImGui::SetNextWindowPos({500.0f, 50.0f}, ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowBgAlpha(0.5f);
    ImGui::Begin("Options", nullptr);
    ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.6f);

    /* General information */
    ImGui::Text("Hide/show menu: H");
    ImGui::Text("Rendering: %3.2f FPS (1 thread)", Double(ImGui::GetIO().Framerate));
    ImGui::Spacing();

    /* Rendering parameters */
    if(ImGui::TreeNodeEx("Particle Rendering", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::PushID("Particle Rendering");
        {
            constexpr const char* items[] = {"Uniform", "Ramp by ID", "Random"};
            static Int colorMode = 1;
            if(ImGui::Combo("Color Mode", &colorMode, items, 3));

            if(colorMode == 0) { /* Uniform color */
                static Color3 color = Color3::blue();
                if(ImGui::ColorEdit3("Diffuse Color", color.data()));
            }
        }
        static Vector3 lightDir = Vector3::xAxis();
        if(ImGui::InputFloat3("Light Direction", lightDir.data()));
        ImGui::PopID();
        ImGui::TreePop();
    }
    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    /* Simulation parameters */
    if(ImGui::TreeNodeEx("Simulation", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::PushID("Simulation");
        Float a,b,c;
        bool d;
        ImGui::InputFloat("Stiffness", &a);
        ImGui::SliderFloat("Viscosity", &b, 0.0f, 1.0f);
        ImGui::SliderFloat("Restitution", &c, 0.0f, 1.0f);
        ImGui::Checkbox("Dynamic Boundary", &d);
        ImGui::PopID();
        ImGui::TreePop();
    }
    ImGui::Spacing();
    ImGui::Separator();

    /* Reset */
    ImGui::Spacing();
    bool p = false;
    if(ImGui::Button(p ? "Play Sim" : "Pause Sim"))
        p ^= true;
    ImGui::SameLine();
    if(ImGui::Button("Reset Sim")) {
        p = false;
    }
    ImGui::SameLine();
    if(ImGui::Button("Reset Camera"));
    ImGui::PopItemWidth();
    ImGui::End();
}

//void Viewer::keyPressEvent(KeyEvent& event) {
//    if(m_imgui.handleKeyPressEvent(event)) return;
//}
//
//void Viewer::keyReleaseEvent(KeyEvent& event) {
//    if(m_imgui.handleKeyReleaseEvent(event)) return;
//}
//
//void Viewer::mousePressEvent(MouseEvent& event) {
//    if(m_imgui.handleMousePressEvent(event)) return;
//}
//
//void Viewer::mouseReleaseEvent(MouseEvent& event) {
//    if(m_imgui.handleMouseReleaseEvent(event)) return;
//}
//
//void Viewer::mouseMoveEvent(MouseMoveEvent& event) {
//    if(m_imgui.handleMouseMoveEvent(event)) return;
//}
//
//void Viewer::mouseScrollEvent(MouseScrollEvent& event) {
//    if(m_imgui.handleMouseScrollEvent(event)) {
//        /* Prevent scrolling the page */
//        event.setAccepted();
//        return;
//    }
//}
//
//void Viewer::textInputEvent(TextInputEvent& event) {
//    if(m_imgui.handleTextInputEvent(event)) return;
//}


//Viewer::Viewer(int argc, char** argv): Platform::Application{Arguments{argc,argv},
//                                                                              Configuration{}.setTitle("Magnum ImGui Example")
//                                                                                      .setWindowFlags(Configuration::WindowFlag::Resizable)}
//{
//    m_imgui = ImGuiIntegration::Context(Vector2{windowSize()}/dpiScaling(),
//                                       windowSize(), framebufferSize());
//
//    /* Set up proper blending to be used by ImGui. There's a great chance
//       you'll need this exact behavior for the rest of your scene. If not, set
//       this only for the drawFrame() call. */
//    GL::Renderer::setBlendEquation(GL::Renderer::BlendEquation::Add,
//                                   GL::Renderer::BlendEquation::Add);
//    GL::Renderer::setBlendFunction(GL::Renderer::BlendFunction::SourceAlpha,
//                                   GL::Renderer::BlendFunction::OneMinusSourceAlpha);
//
//#if !defined(MAGNUM_TARGET_WEBGL) && !defined(CORRADE_TARGET_ANDROID)
//    /* Have some sane speed, please */
//    setMinimalLoopPeriod(16);
//#endif
//}

//void Viewer::drawEvent() {
//    GL::defaultFramebuffer.clear(GL::FramebufferClear::Color);
//
//    bool camChanged = m_camera->update();
//    m_camera->draw(m_scene->drawables());
//
//    m_imgui.newFrame();
//
//    /* Enable text input, if needed */
//    if(ImGui::GetIO().WantTextInput && !isTextInputActive())
//        startTextInput();
//    else if(!ImGui::GetIO().WantTextInput && isTextInputActive())
//        stopTextInput();
//
//    /* 1. Show a simple window.
//       Tip: if we don't call ImGui::Begin()/ImGui::End() the widgets appear in
//       a window called "Debug" automatically */
//    {
//        ImGui::Text("Hello, world!");
//        ImGui::SliderFloat("Float", &_floatValue, 0.0f, 1.0f);
//        if(ImGui::ColorEdit3("Clear Color", _clearColor.data()))
//            GL::Renderer::setClearColor(_clearColor);
//        if(ImGui::Button("Test Window"))
//            _showDemoWindow ^= true;
//        if(ImGui::Button("Another Window"))
//            _showAnotherWindow ^= true;
//        ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
//                    1000.0/Double(ImGui::GetIO().Framerate), Double(ImGui::GetIO().Framerate));
//    }
//
//    /* 2. Show another simple window, now using an explicit Begin/End pair */
//    if(_showAnotherWindow) {
//        ImGui::SetNextWindowSize(ImVec2(500, 100), ImGuiCond_FirstUseEver);
//        ImGui::Begin("Another Window", &_showAnotherWindow);
//        ImGui::Text("Hello");
//        ImGui::End();
//    }
//
//    /* 3. Show the ImGui demo window. Most of the sample code is in
//       ImGui::ShowDemoWindow() */
//    if(_showDemoWindow) {
//        ImGui::SetNextWindowPos(ImVec2(650, 20), ImGuiCond_FirstUseEver);
//        ImGui::ShowDemoWindow();
//    }
//
//    /* Update application cursor */
//    m_imgui.updateApplicationCursor(*this);
//
//    /* Set appropriate states. If you only draw ImGui, it is sufficient to
//       just enable blending and scissor test in the constructor. */
//    GL::Renderer::enable(GL::Renderer::Feature::Blending);
//    GL::Renderer::enable(GL::Renderer::Feature::ScissorTest);
//    GL::Renderer::disable(GL::Renderer::Feature::FaceCulling);
//    GL::Renderer::disable(GL::Renderer::Feature::DepthTest);
//
//    m_imgui.drawFrame();
//
//    /* Reset state. Only needed if you want to draw something else with
//       different state after. */
//    GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
//    GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);
//    GL::Renderer::disable(GL::Renderer::Feature::ScissorTest);
//    GL::Renderer::disable(GL::Renderer::Feature::Blending);
//
//    swapBuffers();
//    redraw();
//}

//void Viewer::viewportEvent(ViewportEvent& event) {
//    GL::defaultFramebuffer.setViewport({{}, event.framebufferSize()});
//
//    m_imgui.relayout(Vector2{event.windowSize()}/event.dpiScaling(),
//                    event.windowSize(), event.framebufferSize());
//}