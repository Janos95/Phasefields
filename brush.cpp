//
// Created by janos on 02.04.20.
//

#include "brush.hpp"
#include "bfs.hpp"
#include "toggle_button.hpp"
#include <scoped_timer/scoped_timer.hpp>

#include <imgui.h>

#include <Magnum/ImageView.h>
#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/PixelFormat.h>
#include <Magnum/Image.h>
#include <Magnum/Math/FunctionsBatch.h>

#include <fmt/core.h>
#include <thread>

using namespace Magnum;
using namespace Corrade;

using namespace std::chrono_literals;


Brush::Brush(PhasefieldData &data): m_phasefieldData(data) {}

void Brush::startPainting() {

    m_thread = std::thread([&]{
        auto& pd = m_phasefieldData;
        auto& vertices = pd.V;
        auto adjacencyList = TriangleMeshAdjacencyList(pd.V, pd.F);
        int step = 1;
        bool reload = false;
        Vector3* it{};
        Vector3 p;
        while(!m_stop){
            if(m_loadPoint) {
                {
                    m_loadPoint = false;
                    std::lock_guard l(*m_mutex);
                    p = m_point;
                }
                it = std::min_element(vertices.begin(), vertices.end(),
                                           [&](auto const &v1, auto const &v2) {
                                               return (v1 - p).dot() < (v2 - p).dot();
                                           });//@todo kdtree here would speed things up considerably
            }
            BreadthFirstSearch bfs(adjacencyList, std::distance(vertices.begin(), it));

            {
                std::lock_guard l(pd.mutex);
                for (int j = 0; (j < step) && bfs.step(); ++j) {
                    for(int a : bfs.enqueuedVertices()){
                        auto& phasefield = pd.phasefield;
                        phasefield[a] = .99 * phasefield[a] + .01f * m_phase;
                    }
                }

                auto textureCoords = pd.meshData.mutableAttribute<Vector2>(Trade::MeshAttribute::TextureCoordinates);

                for (int m = 0; m < textureCoords.size(); ++m) {
                    textureCoords[m].x() = (1.+pd.phasefield[m]) * .5f;
                }
                pd.status = PhasefieldData::Status::PhasefieldUpdated;
            }

            fmt::print("Sleeping for {} seconds\n", (0.01s).count());
            std::this_thread::sleep_for(.01s);
            ++step;
        }
    });
}


Vector3 unproject(Vector2i const& windowPosition, Viewer& viewer){

    auto framebufferSize = viewer.framebufferSize();
    auto windowSize = viewer.windowSize();
    auto const& camera = viewer.camera();

    const Vector2i position = windowPosition*Vector2{framebufferSize}/Vector2{windowSize};
    const Vector2i fbPosition{position.x(), GL::defaultFramebuffer.viewport().sizeY() - position.y() - 1};

    Float depth;
    {
        ScopedTimer timer("Reading depth from framebuffer", true);
        Image2D data = GL::defaultFramebuffer.read(
                Range2Di::fromSize(fbPosition, Vector2i{1}).padded(Vector2i{2}),
                {GL::PixelFormat::DepthComponent, GL::PixelType::Float});
        depth = Math::min<Float>(Containers::arrayCast<const Float>(data.data()));
    }

    const Vector2i viewSize = windowSize;
    const Vector2i viewPosition{windowPosition.x(), viewSize.y() - windowPosition.y() - 1};
    const Vector3 in{2 * Vector2{viewPosition} / Vector2{viewSize} - Vector2{1.0f}, depth * 2.0f - 1.0f};

    //get global coordinates
    return (camera.transformationMatrix() * camera.projectionMatrix().inverted()).transformPoint(in);
}

void Brush::mousePressEvent(Viewer::MouseEvent &event, Viewer& viewer) {
    if(!m_brushing) return;

    m_point = unproject(event.position(), viewer);
    m_stop = false;
    m_loadPoint = true;
    startPainting();

    event.setAccepted();
}

void Brush::mouseMoveEvent(Viewer::MouseMoveEvent& event, Viewer& viewer) {
    if(!m_brushing) return;

    if(m_stop) return;
    auto position = event.position();
    m_loadPoint = true;
    std::lock_guard l(*m_mutex);
    m_point = unproject(position, viewer);
    event.setAccepted();
}

void Brush::mouseReleaseEvent(Viewer::MouseEvent &event, Viewer& viewer) {
    if(!m_brushing) return;

    m_stop = true;
    m_thread.join();
    event.setAccepted();
}

void Brush::drawImGui() {
    if (ImGui::TreeNode("Brush"))
    {
        constexpr int step = 1;
        ImGui::InputScalar("Speed (vertices per second)", ImGuiDataType_U32, &m_speed, &step, nullptr, "%d");

        constexpr float lower = -1.f, upper = 1.f;
        ImGui::SliderScalar("Phase", ImGuiDataType_Float, &m_phase, &lower, &upper, "%.3f", 1.0f);
        if(toggleButton("Enable Brushing", &m_brushing)){
            if(m_brushing)
                GL::defaultFramebuffer.mapForRead(GL::DefaultFramebuffer::ReadAttachment::Front);
            else
                GL::defaultFramebuffer.mapForRead(GL::DefaultFramebuffer::ReadAttachment::None);
        }

        ImGui::TreePop();
    }
}

void Brush::tickEvent(Scene &) {
    auto& pd = m_phasefieldData;
    switch(pd.status){
        case PhasefieldData::Status::NewMesh :
        case PhasefieldData::Status::Subdivided :
            //m_adjacencyList = TriangleMeshAdjacencyList(pd.V, pd.F);
            break;
    }
}