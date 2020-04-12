//
// Created by janos on 02.04.20.
//

#include "brush.hpp"
#include "dijkstra.hpp"
#include "toggle_button.hpp"

#include "geodesic_algorithm_exact.h"

#include <scoped_timer/scoped_timer.hpp>

#include <imgui.h>

#include <Magnum/ImageView.h>
#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/PixelFormat.h>
#include <Magnum/Image.h>
#include <Magnum/Math/FunctionsBatch.h>

#include <fmt/core.h>
#include <thread>
#include <unordered_set>

using namespace Magnum;
using namespace Corrade;

using namespace std::chrono_literals;


Brush::Brush(PhasefieldData &data): m_phasefieldData(data) {}

//not beautiful but its doing the job :)
void Brush::startPainting() {
    using my_dur_t = std::chrono::duration<float, std::ratio<1>>;

    m_thread = std::thread([&]{
        float targetDist = 0.f;
        auto& pd = m_phasefieldData;
        auto& vertices = pd.V;
        auto& faces = pd.F;
        auto& phasefield = pd.phasefield;
        auto recFilter = m_recursiveFilterFactor;
        my_dur_t dur(m_timeTilNextWavefront);
        auto begin = std::chrono::steady_clock::now();
        std::vector<int> history; //@todo use some more efficient data structure here for fast lookup
        bool propagating = true;
        Vector3* it{};
        Vector3 p;

        auto adjacencyList = TriangleMeshAdjacencyList(vertices, faces);
        Dijkstra dijkstra(adjacencyList);

        geodesic::Mesh mesh;
        mesh.initialize_mesh_data(Containers::arrayCast<Float>(vertices), faces);		//create internal mesh data structure including edges
        geodesic::GeodesicAlgorithmExact exactGeodesics1(&mesh), exactGeodesics2(&mesh);
        std::vector<geodesic::SurfacePoint> sources;

        std::vector<std::pair<double, int>> distances(vertices.size());
        while(true){

            if(m_loadPoint.exchange(false)) {
                {
                    std::lock_guard l(m_mutex);
                    p = m_point;
                }
                it = std::min_element(vertices.begin(), vertices.end(),
                                           [&](auto const &v1, auto const &v2) {
                                               return (v1 - p).dot() < (v2 - p).dot();
                                           });//@todo use a kdtree here
                auto source = std::distance(vertices.begin(), it);
                //reset everyhting
                //history.clear();

                //dijkstra.reset(); //@todo this does not need to be done in first iteration...
                //dijkstra.setSource(source);

                ScopedTimer timer("Computing geodesics", true);

                sources.clear();
                sources.emplace_back(&mesh.vertices()[source]);
                exactGeodesics1.propagate(sources);
                for(unsigned i=0; i<mesh.vertices().size(); ++i) {
                    geodesic::SurfacePoint p(&mesh.vertices()[i]);
                    double distance;
                    unsigned best_source = exactGeodesics1.best_source(p, distance);        //for a given surface point, find closets source and distance to this source
                    distances[i].first = distance;
                    distances[i].second = i;
                }
                std::sort(distances.begin(), distances.end(),[](auto& a, auto& b){ return a.first < b.first; });
                    //exactGeodesics2.clear();
                //exactGeodesics2.initialize(sources, geodesic::GEODESIC_INF);
            }

            targetDist += m_distStep;

            for(auto [d, i] : distances) {
                if(d > targetDist) break;
                phasefield[i] = (1.f - recFilter) * phasefield[i] + recFilter * m_phase;
            }

            //for(auto node : history)

            //if(propagating) {
            //    targetDist += m_distStep;
            //    fmt::print("{}\n", targetDist);
            //    while(true){
            //        std::array<double, 2> distances;
            //        std::array<int, 2> nodes{-1,-1};
            //        auto status = exactGeodesics2.step(nodes, distances);
            //        if(geodesic::Status::Continue == status) {
            //            continue;
            //        } else if(geodesic::Status::Done == status){
            //            propagating = false;
            //            break;
            //        } else if(geodesic::Status::NewNode == status){
            //            double d = -1.f;
            //            for (int j = 0; j < 2; ++j) {
            //                auto node = nodes[j];
            //                if(node < 0) continue;
            //                auto dist = distances[j];
            //                {
            //                    double distance;
            //                    geodesic::SurfacePoint p(&mesh.vertices()[node]);
            //                    exactGeodesics1.best_source(p,distance);
            //                    maxDifference = std::max(maxDifference, std::abs(distance - dist));
            //                }

            //                d = std::max(d, dist);
            //                history.push_back(node);
            //                phasefield[node] = (1.f - recFilter) * phasefield[node] + recFilter * m_phase;
            //            }
            //            if(d > targetDist) break;
            //        } else {
            //            CORRADE_INTERNAL_ASSERT(false);
            //        }
            //    }
            //}

            {
                std::lock_guard l(pd.mutex);
                auto textureCoords = pd.meshData.mutableAttribute<Vector2>(Trade::MeshAttribute::TextureCoordinates);

                for (int m = 0; m < textureCoords.size(); ++m) {
                    textureCoords[m].x() = (1.+pd.phasefield[m]) * .5f;
                }
                pd.status = PhasefieldData::Status::PhasefieldUpdated;
            }

            //@todo subtract computation time?
            {
                std::unique_lock ul(m_mutex);
                if(m_stop) return;
                m_cv.wait_for(ul, dur, [this]{ return m_stop; });
                if(m_stop) return;
            }
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
    auto p = unproject(position, viewer);
    std::lock_guard l(m_mutex);
    m_point = p;
    m_loadPoint = true;
    event.setAccepted();
}

void Brush::mouseReleaseEvent(Viewer::MouseEvent &event, Viewer& viewer) {
    if(!m_brushing) return;

    {
        std::lock_guard l(m_mutex);
        m_stop = true;
        m_cv.notify_one();
    }

    m_thread.join();
    event.setAccepted();
}

void Brush::drawImGui() {
    if (ImGui::TreeNode("Brush"))
    {
        constexpr int step = 1;
        constexpr float min = 0.f, max = 1.f;
        ImGui::SliderScalar("Time til next wavefront (in s)", ImGuiDataType_Float, &m_timeTilNextWavefront, &min, &max, "%.3f", 2.0f);
        ImGui::SliderScalar("Recursive Phase Filter Factor", ImGuiDataType_Float, &m_recursiveFilterFactor, &min, &max, "%.3f", 2.0f);
        ImGui::SliderScalar("Distance Step", ImGuiDataType_Float, &m_distStep, &min, &max, "%.5f", 2.0f);
        constexpr float lower = -1.f, upper = 1.f;
        ImGui::SliderScalar("Phase", ImGuiDataType_Float, &m_phase, &lower, &upper, "%.2f", 1.0f);
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