//
// Created by janos on 24.04.20.
//

#include "polygonize_expression.hpp"
#include "primitive_options.hpp"
#include "exprtk.hpp"

#include <Magnum/Math/Vector3.h>
#include <Magnum/Magnum.h>
#include <Corrade/Containers/GrowableArray.h>
#include <Corrade/Utility/Algorithms.h>

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/Surface_mesh.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include <algorithm>

using namespace Magnum;
using namespace Corrade;



// default triangulation for Surface_mesher
using Tr = CGAL::Surface_mesh_default_triangulation_3;
// c2t3
typedef CGAL::Complex_2_in_triangulation_3 <Tr> C2t3;
using GT = Tr::Geom_traits;
using Sphere_3 = GT::Sphere_3;
using Point_3 = GT::Point_3;
using FT = GT::FT;

typedef FT (* Function)(Point_3);

typedef CGAL::Implicit_surface_3 <GT, Function> Surface_3;
using Mesh = CGAL::Surface_mesh<Point_3>;

FT sphere_function(Point_3 p) {
    const FT x2 = p.x()*p.x(), y2 = p.y()*p.y(), z2 = p.z()*p.z();
    return x2 + y2 + z2 - 1;
}


struct ImplicitFunctionFromExpression {
    using symbol_table_t = exprtk::symbol_table<FT>;
    using expression_t = exprtk::expression<FT>;
    using parser_t = exprtk::parser<FT>;

    symbol_table_t symbol_table;
    expression_t expression;
    parser_t parser;

    FT X, Y, Z;

    explicit ImplicitFunctionFromExpression(std::string const& expressionString) {
        symbol_table.add_variable("x", X);
        symbol_table.add_variable("y", Y);
        symbol_table.add_variable("z", Z);
        expression.register_symbol_table(symbol_table);
        parser.compile(expressionString, expression);
    }

    FT operator()(Point_3 p) {
        X = p.x();
        Y = p.y();
        Z = p.z();
        return expression.value();
    }
};

Mg::Trade::MeshData polygonizeExpression(std::string const& expression, PolygonizationOptions const& options) {
    auto expr = expression;
    std::replace(expr.begin(), expr.end(), '\n', ' ');
    Tr tr;            // 3D-Delaunay triangulation
    C2t3 c2t3(tr);   // 2D-complex in 3D-Delaunay triangulation
    // defining the surface
    ImplicitFunctionFromExpression implicitFunction(expr);
    Surface_3 surface([&](Point_3 p) { return implicitFunction(p); },             // pointer to function
                      Sphere_3(CGAL::ORIGIN, options.boundingSphereRadius)); // bounding sphere
    // Note that "2." above is the *squared* radius of the bounding sphere!
    // defining meshing criteria
    CGAL::Surface_mesh_default_criteria_3 <Tr> criteria(options.angleBound,  // angular bound
                                                        options.radiusBound,  // radius bound
                                                        options.distanceBound); // distance bound
    // meshing surface
    CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());

    std::cout << "Final number of points: " << tr.number_of_vertices() << "\n";


    Mesh mesh;
    CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, mesh);

    Containers::Array<Vector3> vertices;
    Containers::Array<UnsignedInt> indices;

    int idx = 0;
    for(auto const& vd : mesh.vertices()){
        auto const& p = mesh.point(vd);
        Containers::arrayAppend(vertices, Containers::InPlaceInit, p.x(), p.y(), p.z());
    }

    for(auto const& fd : mesh.faces()){
        for(auto const& vd : vertices_around_face(mesh.halfedge(fd), mesh)){
            Containers::arrayAppend(indices, vd.idx());
        }
    }

    Trade::MeshIndexData indexData{indices};
    auto indexBuffer = Containers::arrayAllocatorCast<char>(std::move(indices));

    Trade::MeshAttributeData positionData{Trade::MeshAttribute::Position, Containers::arrayCast<Vector3>(vertices)};
    auto positionsBuffer = Containers::arrayAllocatorCast<char>(std::move(vertices));

    return Trade::MeshData{MeshPrimitive::Triangles,
                           std::move(indexBuffer), indexData,
                           std::move(positionsBuffer), {positionData}};
}
