//
// Created by janos on 03.03.20.
//

#include "write_mesh.hpp"

#include "../colormaps.hpp"
#include "../mesh/mesh.hpp"

#include <OpenMesh/Core/IO/MeshIO.hh>

#include <Eigen/Core>

#include <Corrade/Containers/Array.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Magnum.h>

#include <fmt/core.h>

#include <algorithm>

using namespace Corrade;
using namespace Magnum;


auto convertToOpenMesh(
        Containers::ArrayView<const Vector3d> points,
        Containers::ArrayView<const Color4> colors,
        Containers::ArrayView<const Int> faces)
{
    MagnumMesh mesh;
    mesh.request_vertex_colors();

    CORRADE_INTERNAL_ASSERT(points.size() == colors.size());
    std::vector<MagnumMesh::VertexHandle> vertexHandles;
    for (std::size_t i = 0; i < points.size(); ++i) {
        Vector3 p(points[i]);
        auto handle = mesh.add_vertex(p);
        mesh.set_color(handle, colors[i]);
        vertexHandles.push_back(std::move(handle));
    }

    for(auto&& f : Containers::arrayCast<const Vector3i>(faces)){
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
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        const Eigen::VectorXd& U,
        folly::FunctionRef<Color4(double)> mapping)
{
    fmt::print("Writing mesh with to {}\n", path);
    Eigen::MatrixXd VT = V.transpose();
    Eigen::MatrixXi FT = F.transpose();
    auto varr = Containers::arrayCast<const Vector3d>(Containers::ArrayView(VT.data(), V.size()));
    Containers::ArrayView farr(FT.data(), F.size());

    Containers::Array<Color4> color(Containers::NoInit, U.size());

    std::transform(U.begin(), U.end(), color.begin(), mapping);
    auto mesh = convertToOpenMesh(varr, color, farr);

    OpenMesh::IO::write_mesh(mesh, path, OpenMesh::IO::Options::VertexColor | OpenMesh::IO::Options::Binary);
}

