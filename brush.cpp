//
// Created by janos on 02.04.20.
//

#include "brush.hpp"
#include "bfs.hpp"


#include <imgui.h>

#include <Magnum/ImageView.h>
#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/PixelFormat.h>

#include <thread>

using namespace Magnum;
using namespace Corrade;

using namespace std::chrono_literals;

void Brush::paint(UnsignedInt vertexIndex, duration_type dur) {

}

//@todo this is not ready yet
Vector3d unproject(
        Vector3 const& window,
        Matrix4 const& transformation,
        Matrix4 const& projection,
        Vector2i const& viewport){
    Vector4 a {
        2.f * window.x() / viewport.x() - 1.f,
        2.f * window.y() / viewport.y() - 1.f,
        2.f * window.z() - 1.f,
        1.f};
    auto q = (projection * transformation).inverted() * a;
    return Vector3d{q.xyz()};
}

bool Brush::mousePressEvent(Viewer::MouseEvent &event, Viewer& viewer) {
    if(!static_cast<int>(m_phase))
        return false;

    auto& camera = viewer.camera();
    auto position = event.position();
    Float z = -1.f;
    MutableImageView2D depthBufferView{GL::PixelFormat::DepthComponent, GL::PixelType::Float, {1,1}, {&z,1}};
    auto viewport = GL::defaultFramebuffer.viewport();
    GL::defaultFramebuffer.read(viewport, depthBufferView);
    if(z < 0)
        return false;
    auto p = unproject({0,0, z}, camera.transformationMatrix(), camera.projectionMatrix(), camera.viewport());
    auto& vertices = m_phasefieldData->vertices;
    auto it = std::min_element(vertices.begin(), vertices.end(),
            [&](auto const& v1, auto const& v2){ return (v1-p).dot() < (v2-p).dot(); });

    int steps = 1;
    while(m_tracking) {
        BreadthFirstSearch bfs(m_adjacencyList, std::distance(vertices.begin(), it));
        for (int j = 0; j < steps && bfs.step(); ++j) {
            for(int a : bfs.enqueuedVertices()){
                std::lock_guard l(m_phasefieldData->mutex);
                auto& phasefield = m_phasefieldData->phasefield;
                phasefield[a] = .5 * (phasefield[a] + m_phase);
            }
        }
        std::this_thread::sleep_for(1s/60);
        ++steps;
    }
}

bool Brush::mouseMoveEvent(Viewer::MouseMoveEvent &event, Viewer &) {
    if(!m_tracking) return false;
    std::lock_guard l(m_mutex);
    m_position = event.position();
}

bool Brush::mouseReleaseEvent(Viewer::MouseEvent &event, Viewer& viewer) {
    m_tracking = false;
}

void Brush::drawImGui() {
    if (ImGui::TreeNode("Brush"))
    {
        constexpr int step = 1;
        ImGui::InputScalar("Speed (vertices per second)", ImGuiDataType_U32, &m_speed, &step, nullptr, "%d");
        ImGui::InputScalar("Speed (vertices per second)", ImGuiDataType_U32, &m_speed, &step, nullptr, "%d");
        ImGui::TreePop();
    }
}