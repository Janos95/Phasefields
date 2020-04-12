//
// Created by janos on 02.04.20.
//

#pragma once

#include "viewer.hpp"
#include "phasefield_data.hpp"

#include <Eigen/SparseCore>

#include <thread>
#include <condition_variable>

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
        std::vector<Eigen::Triplet<float>> triplets;
        for(auto const& face : faces){
            for (int j = 0; j < 3; ++j) {
                int a = face[j], b = face[(j+1)%3];
                auto w = (vertices[a] - vertices[b]).length();
                triplets.emplace_back(a,b,w);
                triplets.emplace_back(b,a,w);
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

class Brush : public Viewer::AbstractEventHandler {
public:
    using duration_type = std::chrono::duration<double>;

    void startPainting();
    void stopPainting();

    explicit Brush(PhasefieldData& data);

    void mousePressEvent(Viewer::MouseEvent& event, Viewer&) override;
    void mouseMoveEvent(Viewer::MouseMoveEvent& event, Viewer&) override;
    void mouseReleaseEvent(Viewer::MouseEvent& event, Viewer&) override;
    void drawImGui() override;
    void tickEvent(Scene&) override;


private:
    PhasefieldData& m_phasefieldData;

    Magnum::UnsignedInt m_speed = 1;
    float m_phase = 0;
    float m_recursiveFilterFactor = 0.1;
    float m_timeTilNextWavefront = .01;
    float m_distStep = 0.01;
    bool m_brushing = false;


    std::thread m_thread;
    std::condition_variable m_cv;
    bool m_stop = true;

    std::mutex m_mutex;
    std::atomic_bool m_loadPoint = false;
    Vector3 m_point;
};



