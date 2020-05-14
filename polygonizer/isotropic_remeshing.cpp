//
// Created by janos on 12.05.20.
//
#include "isotropic_remeshing.hpp"

#include <Magnum/Math/Vector3.h>
#include <Corrade/Containers/GrowableArray.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <boost/function_output_iterator.hpp>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
using Point = K::Point_3;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor     edge_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;

using namespace Magnum;
using namespace Corrade;

struct halfedge2edge
{
    halfedge2edge(const Mesh& m, std::vector<edge_descriptor>& edges)
            : m_mesh(m), m_edges(edges)
    {}
    void operator()(const halfedge_descriptor& h) const
    {
        m_edges.push_back(edge(h, m_mesh));
    }
    const Mesh& m_mesh;
    std::vector<edge_descriptor>& m_edges;
};

using Index = Mesh::Vertex_index;

void isotropicRemeshing(
        Containers::Array<Vector3d>& vertices,
        Containers::Array<UnsignedInt>& indices)
{
    auto triangles = Containers::arrayCast<Vector3ui>(indices);
    double target_edge_length = 0.01;
    unsigned int nb_iter = 3;

    Mesh mesh;
    std::vector<Index> vs;
    vs.reserve(vertices.size());
    for(auto const& v : vertices){
        vs.push_back(mesh.add_vertex(Point(v.x(), v.y(), v.z())));
    }
    for(auto const& t : triangles){
        mesh.add_face(vs[t[0]], vs[t[1]], vs[t[2]]);
    }

    std::vector<edge_descriptor> border;
    PMP::border_halfedges(faces(mesh), mesh, boost::make_function_output_iterator(halfedge2edge(mesh, border)));
    PMP::split_long_edges(border, target_edge_length, mesh);
    PMP::isotropic_remeshing(
            faces(mesh),
            target_edge_length,
            mesh,
            PMP::parameters::number_of_iterations(nb_iter)
                    .protect_constraints(true)//i.e. protect border, here
    );

    Containers::arrayResize(vertices, Containers::NoInit, mesh.num_vertices());
    Containers::arrayResize(indices, Containers::NoInit, mesh.num_faces() * 3);
    triangles = Containers::arrayCast<Vector3ui>(indices);

    std::transform(mesh.vertices_begin(), mesh.vertices_end(), vertices.begin(), [&](auto& vd) {
        auto const &p = mesh.point(vd);
        return Vector3d(p.x(), p.y(), p.z());
    });

    std::transform(mesh.faces_begin(), mesh.faces_end(), triangles.begin(), [&](auto const& fd) {
        auto incidents = CGAL::vertices_around_face(mesh.halfedge(fd), mesh);
        Vector3ui t;
        std::transform(incidents.begin(), incidents.end(), t.data(), [](auto const& vd){ return vd.idx(); });
        return t;
    });
}