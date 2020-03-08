//
// Created by janos on 02.03.20.
//


#pragma once

#include "mesh.hpp"
#include <Corrade/Containers/ArrayView.h>

using namespace Magnum;
using namespace Corrade;

auto convertToOpenMesh(
        Containers::ArrayView<Vector3> points,
        Containers::ArrayView<Color4> colors,
        Containers::ArrayView<UnsignedInt> faces)
{
    MagnumMesh mesh;
    mesh.request_vertex_colors();

    CORRADE_INTERNAL_ASSERT(points.size() == colors.size());
    std::vector<MagnumMesh::VertexHandle> vertexHandles;
    for (int i = 0; i < points.size(); ++i) {
       auto handle = mesh.add_vertex(points[i]);
       mesh.set_color(handle, colors[i]);
       vertexHandles.push_back(std::move(handle));
    }

    for(auto&& f : Containers::arrayCast<Vector3ui>(faces)){
        std::vector<MagnumMesh::VertexHandle> face;
        for (int i = 0; i < 3; ++i) {
            face.push_back(vertexHandles[f[i]]);
        }
        mesh.add_face(face);
    }

   return mesh;
}
