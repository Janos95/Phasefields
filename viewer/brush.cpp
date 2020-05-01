//
// Created by janos on 02.04.20.
//

#include "brush.hpp"
#include "dijkstra.hpp"
#include "upload.hpp"

#include <scoped_timer/scoped_timer.hpp>

#include <imgui.h>

#include <Magnum/ImageView.h>
#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/PixelFormat.h>
#include <Magnum/Image.h>
#include <Magnum/Math/FunctionsBatch.h>

#include <Eigen/SparseCore>

#include <Corrade/Utility/Algorithms.h>

using namespace Magnum;
using namespace Corrade;

using namespace std::chrono_literals;


enum class Side : UnsignedByte {
    Left,
    Right
};

template<Side side, UnsignedInt Size, class T>
Math::Range<Size, T> trim(
        Math::Range<Size, T> range,
        UnsignedInt dim,
        typename Math::Range<Size,T>::VectorType const& p)
{
    if constexpr(side == Side::Right) range.max()[dim] = p[dim];
    else range.min()[dim] = p[dim];
    return range;
}

template<class Vector>
auto distanceSq(Vector const& p, Math::Range<Vector::Size, typename Vector::Type> const& range){
    typename Vector::Type dsq{0};
    for (int j = 0; j < Vector::Size; ++j){
        auto lower = range.min()[j];
        auto upper = range.max()[j];
        if(p[j] < lower)
            dsq += std::pow(lower - p[j], 2);
        else if(p[j] > upper)
            dsq += std::pow(p[j] - upper, 2);
    }
    return dsq;
}

template<class Vector>
class KDTree {
public:

    struct Node {
        constexpr static auto Invalid = -1;
        int leftChild;
        int rightChild;
        int pointIndex;
    };

    using Type = typename Vector::Type;
    enum {
        Size = Vector::Size
    };
    using Range = Math::Range<Size, Type>;

    explicit KDTree(Containers::ArrayView<Vector> const& points) :
        m_nodes(Containers::NoInit, points.size()),
        m_points(points)
    {
        Containers::Array<UnsignedInt> nodeIds(Containers::NoInit, points.size());
        std::iota(nodeIds.begin(), nodeIds.end(), 0);
        int size = 0;
        construct(nodeIds, 0, size);
    }

    Vector findMin(int axis, int cd){

    }

    struct NNResult{
        int pointIndex;
        Type distanceSquared = std::numeric_limits<Type>::infinity();
    };

    /**
     * @param queryPoint the point for which to perform
     * the nearest neighbor query.
     * @return struct containing index corresponding to nearst neighbor and
     * squared distance
     */
    NNResult nearestNeighbor(Vector const& queryPoint){
        NNResult result;
        recurse(queryPoint, m_root, 0, m_bb, result);
        return result;
    }

private:

    int construct(Containers::ArrayView<UnsignedInt> const& ids, int depth, int& size){
        if(ids.empty()) return Node::Invalid;

        auto b = ids.begin(), e = ids.end();
        UnsignedInt sizeFront = (e-b)/2;
        UnsignedInt sizeBack = ids.size() - sizeFront;
        auto m = b + sizeFront;
        auto cd = depth % Size;
        std::nth_element(b, m, e, [&](auto id1, auto id2){ return m_points[id1][cd] < m_points[id2][cd]; });

        auto nodeIdx = size;
        auto& node = m_nodes[size++];
        node.pointIndex = *m;
        node.leftChild = construct({b, sizeFront}, depth + 1, size);
        node.rightChild = construct({m, sizeBack}, depth + 1, size);
        return nodeIdx;
    }

    void recurse(Vector const& q, int nodeId, int cd, Range const& bb, NNResult& result){
        auto const& node = m_nodes[nodeId];
        auto const& p = m_points[node.pointIndex];
        if(nodeId == Node::Invalid || distanceSq(p, bb) > result.distanceSquared) return;
        auto distSq = (p - q).dot();
        if(distSq < result.distanceSquared){
            result.distanceSquared = distSq;
            result.pointIndex = node.pointIndex;
        }
        auto nextCd = (cd+1) % Size;
        if(q[cd] < p[cd]){
            recurse(q, node.leftChild, nextCd, trim<Side::Left>(bb, cd, p), result);
            recurse(q, node.rightChild, nextCd, trim<Side::Right>(bb, cd, p), result);
        } else {
            recurse(q, node.rightChild, nextCd, trim<Side::Right>(bb, cd, p), result);
            recurse(q, node.leftChild, nextCd, trim<Side::Left>(bb, cd, p), result);
        }
    }

    int m_root;
    Range m_bb;
    Containers::Array<Node> m_nodes;
    Containers::ArrayView<Vector> m_points;
};


struct TriangleMeshAdjacencyList
{
    using graph_type = Eigen::SparseMatrix<float>;
    using sparse_iterator_type = typename graph_type::InnerIterator;
    graph_type graph;
    TriangleMeshAdjacencyList() = default;
    TriangleMeshAdjacencyList(Containers::ArrayView<Vector3> vertices, Containers::ArrayView<UnsignedInt> indices):
            graph(vertices.size(), vertices.size())
    {
        auto faces = Containers::arrayCast<Vector3ui>(indices);
        Containers::Array<Eigen::Triplet<float>> triplets;
        for(auto const& face : faces){
            for (int j = 0; j < 3; ++j) {
                int a = face[j], b = face[(j+1)%3];
                auto w = (vertices[a] - vertices[b]).length();
                Containers::arrayAppend(triplets, Containers::InPlaceInit, a, b, w);
                Containers::arrayAppend(triplets, Containers::InPlaceInit, b, a, w);
            }
        }
        graph.setFromTriplets(triplets.begin(), triplets.end(), [](auto const& a, auto const& b){ return b; });
        graph.makeCompressed();
    }

    auto operator[](int vertexId) const{
        return VertexRange{
                Iterator{graph.valuePtr(), graph.innerIndexPtr(), graph.outerIndexPtr()[vertexId]},
                Iterator{graph.valuePtr(), graph.innerIndexPtr(), graph.outerIndexPtr()[vertexId + 1]}
        };
    }

    [[nodiscard]] auto size() const { return graph.rows();}

    struct Iterator {
        Iterator& operator++() {
            ++idx;
            return *this;
        }

        bool operator !=(Iterator const& other) const {
            return idx != other.idx;
        }

        auto operator*() const {
            return std::pair(rows[idx], values[idx]);
        }

        float const* values;
        int const* rows;
        int idx;
    };

    struct VertexRange{
        Iterator m_begin, m_end;
        [[nodiscard]] auto begin() const { return m_begin; }
        [[nodiscard]] auto end() const { return m_end; }
    };
};



Vector3 unproject(Vector2i const& windowPosition, Viewer& viewer){

    auto framebufferSize = viewer.framebufferSize();
    auto windowSize = viewer.windowSize();
    auto const& camera = *viewer.camera;

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

void geodesicSearch(
        Containers::Array<std::pair<Double, int>>& distances,
        Vector3d const& point,
        Containers::Array<Vector3d> const& vertices,
        Containers::Array<UnsignedInt> const& triangles){

    geodesic::Mesh mesh;
    mesh.initialize_mesh_data(Containers::arrayCast<const Double>(vertices), triangles);		//create internal mesh data structure including edges
    geodesic::GeodesicAlgorithmExact exactGeodesics(&mesh);

    Containers::arrayResize(distances, vertices.size());
    auto it = std::min_element(vertices.begin(), vertices.end(),
                                  [&](auto const &v1, auto const &v2) {
                                      return (v1 - point).dot() < (v2 - point).dot();
                                  });
    auto source = it - vertices.begin();
    ScopedTimer timer("Computing geodesics", true);

    std::vector<geodesic::SurfacePoint> sources;
    sources.clear();
    sources.emplace_back(&mesh.vertices()[source]);
    exactGeodesics.propagate(sources);
    for(unsigned i=0; i<mesh.vertices().size(); ++i) {
        geodesic::SurfacePoint p(&mesh.vertices()[i]);
        double distance;
        exactGeodesics.best_source(p, distance);
        distances[i].first = distance;
        distances[i].second = i;
    }
    std::sort(distances.begin(), distances.end(),[](auto& a, auto& b){ return a.first < b.first; });
}

void Brush::mousePressEvent(Viewer::MouseEvent &event, Viewer& viewer) {
    if(!brushing) return;
    point = Vector3d(unproject(event.position(), viewer));
    mesh.initialize_mesh_data(Containers::arrayCast<Double>(pd.V), pd.F);		//create internal mesh data structure including edges
    geodesicSearch();
    stop = false;
    event.setAccepted();
}

void Brush::mouseMoveEvent(Viewer::MouseMoveEvent& event, Viewer& viewer) {
    if(!brushing || stop) return;
    auto position = event.position();
    point = Vector3d(unproject(position, viewer));
    geodesicSearch();
    event.setAccepted();
}

void Brush::mouseReleaseEvent(Viewer::MouseEvent &event, Viewer& viewer) {
    stop = true;
    if(!brushing) return;
    event.setAccepted();
}

void Brush::keyPressEvent(Viewer::KeyEvent& event, Viewer&){
    if(event.key() == Viewer::KeyEvent::Key::LeftShift){
        brushing = true;
        GL::defaultFramebuffer.mapForRead(GL::DefaultFramebuffer::ReadAttachment::Front);
    }
}
void Brush::keyReleaseEvent(Viewer::KeyEvent& event, Viewer&){
    if(event.key() == Viewer::KeyEvent::Key::LeftShift){
        brushing = false;
        GL::defaultFramebuffer.mapForRead(GL::DefaultFramebuffer::ReadAttachment::None);
    }
}

void Brush::drawImGui(Viewer&) {
    if (ImGui::TreeNode("Brush"))
    {
        constexpr int step = 1;
        constexpr double stepDist = 0.01;
        constexpr double min = 0.f, max = 1.f;
        ImGui::SliderScalar("Recursive Phase Filter Factor", ImGuiDataType_Double, &recursiveFilterFactor, &min, &max, "%.3f", 2.0f);
        ImGui::SliderScalar("Distance Step", ImGuiDataType_Double, &distStep, &min, &max, "%.5f", 2.0f);
        ImGui::InputScalar("Maximal Distance", ImGuiDataType_Double, &maxDist, &stepDist, nullptr, "%.3f");
        constexpr double lower = -1.f, upper = 1.f;
        ImGui::SliderScalar("Phase", ImGuiDataType_Double, &phase, &lower, &upper, "%.2f", 1.0f);

        const auto colBrushing = ImVec4(0.56f, 0.83f, 0.26f, 1.0f);
        const auto colNotBrushing = ImVec4(0.85f, 0.85f, 0.85f, 1.0f);
        ImGui::TextColored(brushing ? colBrushing : colNotBrushing, "Press (Left) Shift To Enable Brushing");
        ImGui::SameLine();
        ImVec2 p = ImGui::GetCursorScreenPos();
        ImDrawList* draw_list = ImGui::GetWindowDrawList();
        float height = ImGui::GetFrameHeight();
        float width = height * 1.55f;
        ImGui::InvisibleButton("brush modus", ImVec2(width, height));
        draw_list->AddRectFilled(p, ImVec2(p.x + width, p.y + height), brushing ? ImGui::GetColorU32(colBrushing) : ImGui::GetColorU32(colNotBrushing), height * 0.1f);
        ImGui::TreePop();
    }
}
