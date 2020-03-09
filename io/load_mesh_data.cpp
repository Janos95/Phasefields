//
// Created by janos on 2/8/20.
//
#include <scoped_timer/scoped_timer.hpp>
#include "load_mesh_data.hpp"

#include <Corrade/Containers/Pointer.h>
#include <Corrade/Containers/Optional.h>

#include <MagnumPlugins/AssimpImporter/AssimpImporter.h>

using namespace Magnum;
using namespace Corrade;

Trade::MeshData loadMeshData(const std::string& path){

    PluginManager::Manager<Trade::AbstractImporter> manager;
    auto importer = manager.loadAndInstantiate("AssimpImporter");
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

