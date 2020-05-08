//
// Created by janos on 08.05.20.
//
#include "polygonize_wireframe.hpp"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_segment_primitive.h>

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>

#include <Corrade/Containers/GrowableArray.h>
#include <Magnum/Math/Vector3.h>
#include <Magnum/Magnum.h>

using namespace Magnum;
using namespace Corrade;

using Scalar = double;
using K = CGAL::Simple_cartesian<Scalar>;
using Point = K::Point_3;
using Sphere = K::Sphere_3;

using Triangulation = CGAL::Surface_mesher::Surface_mesh_default_triangulation_3_generator<K>::Type;
using C2t3 = CGAL::Complex_2_in_triangulation_3<Triangulation>;
using Traits = Triangulation::Geom_traits;
using ImplicitSurface = CGAL::Implicit_surface_3<K>;
using Mesh = CGAL::Surface_mesh<Point>;

using Segment = K::Segment_3;
using Iterator = std::vector<Segment>::iterator;
using Primitive = CGAL::AABB_segment_primitive<K, Iterator>;
using AABBTraits = CGAL::AABB_traits<K, Primitive>;
using Tree = CGAL::AABB_tree<AABBTraits>;

using Mesh = CGAL::Surface_mesh<Point>;

using namespace CGAL::parameters;

Mg::Trade::MeshData polygonizeWireframe(
        Cr::Containers::ArrayView<const Mg::Vector3d> const& vertices,
        Cr::Containers::ArrayView<const Mg::Vector3ui> const& triangles,
        PolygonizationOptions const& options)
{
    std::vector<std::pair<UnsignedInt, UnsignedInt>> edges;
    edges.reserve(3 * triangles.size());

    for (int i = 0; i < triangles.size(); ++i) {
        for (int j = 0; j < 3; ++j){
            edges.emplace_back(triangles[i][j], triangles[i][(j+1)%3]);
        }
    }

    std::sort(edges.begin(), edges.end());
    auto it = std::unique(edges.begin(), edges.end());
    edges.erase(it, edges.end());

    std::vector<Segment> segments(edges.size());
    for (int i = 0; i < edges.size(); ++i) {
        auto& a = vertices[edges[i].first];
        auto& b = vertices[edges[i].second];
        segments[i] = Segment(Point(a.x(), a.y(), a.z()), Point(b.x(), b.y(), b.z()));
    }

    // constructs the AABB tree and the internal search tree for
    // efficient distance computations.
    Tree tree(segments.begin(),segments.end());
    tree.accelerate_distance_queries();
    // counts #intersections with a plane query

    const Scalar thicknessSq = 0.01 * 0.01;
    auto implicitFunction = [&](Point p){ return tree.squared_distance(p) - thicknessSq; };

    Triangulation tr;            // 3D-Delaunay triangulation
    C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation
    // defining the surface
    ImplicitSurface surface(implicitFunction,             // pointer to function
                      Sphere(CGAL::ORIGIN, options.boundingSphereRadius)); // bounding sphere
    // Note that "2." above is the *squared* radius of the bounding sphere!
    // defining meshing criteria
    CGAL::Surface_mesh_default_criteria_3<Triangulation> criteria(options.angleBound,  // angular bound
                                                       options.radiusBound,  // radius bound
                                                       options.distanceBound); // distance bound
    // meshing surface
    CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());

    std::cout << "Final number of points: " << tr.number_of_vertices() << "\n";


    Mesh mesh;
    CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, mesh);


    Containers::Array<Vector3> vs;
    Containers::Array<UnsignedInt> is;

    int idx = 0;
    for(auto const& vd : mesh.vertices()){
        auto const& p = mesh.point(vd);
        Containers::arrayAppend(vs, Containers::InPlaceInit, p.x(), p.y(), p.z());
    }

    for(auto const& fd : mesh.faces()){
        for(auto const& vd : vertices_around_face(mesh.halfedge(fd), mesh)){
            Containers::arrayAppend(is, vd.idx());
        }
    }

    Trade::MeshIndexData indexData{is};
    auto indexBuffer = Containers::arrayAllocatorCast<char>(std::move(is));

    Trade::MeshAttributeData positionData{Trade::MeshAttribute::Position, Containers::arrayCast<Vector3>(vs)};
    auto positionsBuffer = Containers::arrayAllocatorCast<char>(std::move(vs));

    return Trade::MeshData{MeshPrimitive::Triangles,
                           std::move(indexBuffer), indexData,
                           std::move(positionsBuffer), {positionData}};

}
