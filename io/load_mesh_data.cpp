//
// Created by janos on 2/8/20.
//
#include <scoped_timer/scoped_timer.hpp>
#include "load_mesh_data.hpp"

#include <Eigen/Core>

#include <Corrade/Containers/Pointer.h>
#include <Corrade/Containers/Optional.h>
#include <Corrade/PluginManager/PluginMetadata.h>

#include <Magnum/Math/Vector3.h>
#include <MagnumPlugins/AssimpImporter/AssimpImporter.h>
#include <Magnum/EigenIntegration/Integration.h>

#include <fmt/core.h>

using namespace Magnum;
using namespace Corrade;

Trade::MeshData loadMeshData(const std::string& path){

    PluginManager::Manager<Trade::AbstractImporter> manager;
    auto importer = manager.loadAndInstantiate("AssimpImporter");
    auto name = importer->metadata()->name();
    fmt::print("Trying to load mesh using {}\n", name);
    if(!importer) std::exit(1);

    Debug{} << "Opening file" << path.c_str();

    if(!importer->openFile(path)){
        puts("could not open file");
        std::exit(4);
    }

    Debug{} << "Imported " << importer->meshCount() << " meshes";

    if(!importer->meshCount()){
        puts("No mesh imported, exiting");
        std::exit(5);
    }

    return *importer->mesh(0);
}


std::tuple<Eigen::MatrixXd, Eigen::MatrixXi> toEigen(Trade::MeshData const& mesh){
    auto varr = mesh.positions3DAsArray();
    auto farr = mesh.indicesAsArray();

    using MatrixXU = Eigen::Matrix<UnsignedInt, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    Eigen::MatrixXi F = Eigen::Map<const MatrixXU>(farr.data(), farr.size() / 3, 3).cast<int>();
    Eigen::MatrixXd V(varr.size(), 3);
    for (std::size_t i = 0; i < mesh.vertexCount(); ++i) {
        V.row(i) = EigenIntegration::cast<Eigen::Vector3f>(varr[i]).cast<double>();
    }
    return std::make_tuple(std::move(V), std::move(F));
}
