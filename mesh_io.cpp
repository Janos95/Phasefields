//
// Created by janos on 03.04.20.
//

#include "mesh_io.hpp"

#include <imgui.h>
#include <MagnumPlugins/AssimpImporter/AssimpImporter.h>

#include <assimp/Exporter.hpp>
#include <assimp/scene.h>


void MeshIO::drawEvent() {
    if (ImGui::TreeNode("Mesh IO"))
    {
        ImGui::InputText("Mesh Path (Loading)", m_inputPath.data(), m_inputPath.size());
        ImGui::InputText("Mesh Path (Saving)", m_outputPath.data(), m_outputPath.size());
        if(ImGui::Button("load mesh")){
            loadMesh(std::string{m_inputPath.begin(),std::find(m_inputPath.begin(), m_inputPath.end(), '\0')});
        }
        ImGui::TreePop();
    }
}

template <class To, class From>
To bit_cast(const From &src) noexcept
{
    To dst;
    std::memcpy(&dst, &src, sizeof(To));
    return dst;
}

bool MeshIO::saveMesh(std::string const& path) {
    //@todo generate path if not existent
    //@todo synchronization


    auto* mesh = new aiMesh;
    if(m_meshData->isIndexed()){
        CORRADE_INTERNAL_ASSERT(m_meshData->primitive() == MeshPrimitive::Triangles);
        auto numTriangles = m_meshData->indexCount() / 3;
        mesh->mNumFaces = numTriangles;
        auto faces = new aiFace[numTriangles];
        for (std::size_t i = 0; i < numTriangles; ++i) {
            auto indices = new unsigned int[3];
            for (int j = 0; j < 3; ++j) {
                indices[j] = m_F(i,j);
            }
            faces[i].mIndices = indices;
            faces[i].mNumIndices = 3;
        }
        mesh->mFaces = faces;
    }

    auto numVertices = m_V.rows();
    mesh->mNumVertices = numVertices;

    auto vertexView = m_meshData->attribute<Vector3>(Trade::MeshAttribute::Position);
    auto normalView = m_meshData->attribute<Vector3>(Trade::MeshAttribute::Position);
    auto colorView = m_meshData->attribute<Color4>(Trade::MeshAttribute::Position);
    auto vertices = new aiVector3D[m_V.rows()];
    auto normals = new aiVector3D[m_V.rows()];
    auto colors = new aiColor4D[m_V.rows()];

    std::transform(vertexView.begin(), vertexView.end(), vertices, &bit_cast<aiVector3D, Vector3>);
    std::transform(normalView.begin(), normalView.end(), normals, &bit_cast<aiVector3D, Vector3>);
    std::transform(colorView.begin(), colorView.end(), colors, &bit_cast<aiColor4D, Color4>);

    mesh->mColors[0] = colors;
    mesh->mVertices = vertices;
    mesh->mNormals = normals;

    aiScene scene;
    auto meshes = new aiMesh*;
    meshes[0] = mesh;
    scene.mMeshes = meshes;

    Assimp::Exporter exporter;
    auto desc = exporter.GetExportFormatDescription(0);
    switch (exporter.Export(&scene, desc->id, path.c_str(), 0)) {
        case aiReturn_SUCCESS : return true;
        case aiReturn_OUTOFMEMORY : std::exit(137);
        default : return false;
    }
}


bool MeshIO::loadMesh(std::string const& path){
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

    if(importer->meshCount()){

        if(m_optimizationFuture.valid() && !m_optimizationFuture.isReady()){
            m_exit = true;
            m_optimizationFuture.wait();
        }

        m_meshData = importer->mesh(0);
        if(!m_meshData) return false;
        *m_meshData = preprocess(*m_meshData, CompileFlag::GenerateSmoothNormals);
        std::tie(m_V,m_F) = toEigen(*m_meshData);
        m_U = Eigen::VectorXd::Zero(m_V.rows());
        m_cost->clear();
        m_cost->emplace_back<InterfaceEnergy>("Modica Mortola", 1., m_V, m_F, m_epsilon);
        m_cost->emplace_back<PotentialEnergy>("Modica Mortola", 1., m_V, m_F, m_epsilon);
        m_cost->emplace_back<AreaRegularizer>("Area Regularization", 1., m_V, m_F);
        m_cost->emplace_back<ConnectednessConstraint>("Connectedness Constraint Positive", 1., m_V, m_F, m_epsilon, m_posA, m_posB);
        m_cost->emplace_back<ConnectednessConstraint>("Connectedness Constraint Negative", 1., m_V, m_F, m_epsilon, -m_negA, -m_negB);
        m_problem.emplace(m_cost);
        Debug{} << "Succeded importing mesh";
        m_newMeshData = true;
        return true;
    }
    else{
        Debug{} << "Could not load mesh";
        return false;
    }
}
