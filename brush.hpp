//
// Created by janos on 02.04.20.
//

#pragma once

#include "viewer.hpp"
#include "phasefield_data.hpp"

#include <Eigen/SparseCore>

struct TriangleMeshAdjacencyList
{
    using graph_type = Eigen::SparseMatrix<int>;
    using sparse_iterator_type = typename graph_type::InnerIterator;
    graph_type graph;

    TriangleMeshAdjacencyList(Containers::ArrayView<Vector3d> vertices, Containers::ArrayView<Vector3i> faces):
            graph(vertices.size(), vertices.size())
    {
        std::vector<Eigen::Triplet<int>> triplets;
        for(auto const& face : faces){
            for (int j = 0; j < 3; ++j) {
                int a = face[j], b = face[(j+1)%3];
                triplets.emplace_back(a,b,1);
                triplets.emplace_back(b,a,1);
            }
        }
        graph.setFromTriplets(triplets.begin(), triplets.end(), [](auto const& a, auto const& b){ return b; });
        graph.makeCompressed();
    }

    auto operator[](int vertexId){
        return VertexRange{Iterator{graph, vertexId}, Iterator{graph, vertexId /*@todo */}};
    }

    struct Iterator : public sparse_iterator_type {
        bool operator !=(Iterator const& other){
            //@todo
        }

        auto operator*(){
            return std::pair(this->row(), this->value());
        }
    };

    struct VertexRange{
        Iterator m_begin, m_end;
        auto begin() { return m_begin; }
        auto end() { return m_end; }
    };
};

class Brush : public Viewer::AbstractEventHandler {

    using duration_type = std::chrono::duration<double>;

    void paint(Magnum::UnsignedInt vertexIndex, duration_type dur);

    Brush(PhasefieldData* data):
        m_phasefieldData(data),
        m_adjacencyList(data->vertices, data->faces)
    {
    }

    bool mousePressEvent(Viewer::MouseEvent& event, Viewer&) override;
    bool mouseMoveEvent(Viewer::MouseMoveEvent& event, Viewer&) override;
    bool mouseReleaseEvent(Viewer::MouseEvent& event, Viewer&) override;
    void drawImGui() override;

    PhasefieldData* m_phasefieldData;

    TriangleMeshAdjacencyList m_adjacencyList;
    Magnum::UnsignedInt m_speed;
    float m_phase = 0;
    Vector2i m_position;
    std::atomic_bool m_tracking = false;
    std::mutex m_mutex;
};



