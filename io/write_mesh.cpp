//
// Created by janos on 03.03.20.
//


#include "../colormaps.hpp"
#include "../mesh/mesh.hpp"

#include <OpenMesh/Core/IO/MeshIO.hh>

#include <Eigen/Core>

#include <Corrade/Containers/Array.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Magnum.h>

#include <algorithm>

using namespace Corrade;
using namespace Magnum;


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

void writeMesh(
        const std::string& path,
        Containers::ArrayView<Vector3> varr,
        Containers::ArrayView<UnsignedInt> farr,
        const Eigen::VectorXd& U,
        bool normalize){

    Containers::Array<Color4> color(Containers::NoInit, U.size());
    Eigen::VectorXd V = U;
    if(normalize) {
        double min = U.minCoeff();
        double max = U.maxCoeff();
        printf("(min,max): (%f, %f)\n", min, max);
        V = (U.array() - min) * 1. / (max - min);
    } else{
        V = (V.array() + 1.) / 2.;
    }
    std::transform(V.begin(), V.end(), color.begin(), jet_colormap);
    auto mesh = convertToOpenMesh(varr, color, farr);

    OpenMesh::IO::write_mesh(mesh, path, OpenMesh::IO::Options::VertexColor | OpenMesh::IO::Options::Binary);
}

