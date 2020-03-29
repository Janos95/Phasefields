//
// Created by janos on 26.03.20.
//

#include "optimization_context.hpp"
#include "colormaps.hpp"
#include "modica_mortola.hpp"
#include "connectedness_constraint.hpp"
#include "upload.hpp"

#include <Corrade/Containers/Pointer.h>
#include <Corrade/Containers/Optional.h>
#include <Corrade/PluginManager/PluginMetadata.h>
#include <Corrade/Utility/Algorithms.h>

#include <Magnum/Math/Vector3.h>
#include <MagnumPlugins/AssimpImporter/AssimpImporter.h>
#include <Magnum/EigenIntegration/Integration.h>
#include <Magnum/ImGuiIntegration/Context.hpp>

#include <assimp/Exporter.hpp>
#include <assimp/scene.h>

#include <fmt/core.h>

using namespace Magnum;
using namespace Corrade;


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

OptimizationContext::OptimizationContext():
    m_optimizationCallback(&OptimizationContext::optimizationCallback),
    m_options{.max_num_iterations = 1, .callbacks = {&m_optimizationCallback}},
    m_optThread([this]{ work(); }),
    m_cost(new SumProblem()),
    m_problem(m_cost)
{
}






ceres::CallbackReturnType OptimizationContext::optimizationCallback(ceres::IterationSummary const&){
    //@todo synchronization
    auto colorView = m_meshData->mutableAttribute<Color4>(Trade::MeshAttribute::Color);
    std::transform(m_U.begin(), m_U.end(), colorView.begin(), JetColorMap{});
    m_reupload = true;
    if(!m_run)
        return ceres::SOLVER_TERMINATE_SUCCESSFULLY;
    return ceres::SOLVER_CONTINUE;
}

void OptimizationContext::updateScene(Scene* scene){
    std::lock_guard l(m_meshDataMutex);

    if(m_newMeshData) {
        if(scene) scene->reset();
        scene->addObject("mesh", upload(*m_meshData));
        scene->addNode("mesh", "mesh", ShaderType::Phong);
        scene->setDirty();
    }

    if(m_reupload){
        auto* obj = scene->getObject("mesh");
        CORRADE_INTERNAL_ASSERT(obj);
        GL::Buffer& vertices = obj->vertices;
        auto data = vertices.map(0,
                                 vertices.size(),
                                 GL::Buffer::MapFlag::Write);

        CORRADE_CONSTEXPR_ASSERT(data, "could not map vertex data");
        Utility::copy(m_meshData->vertexData(), data);
        vertices.unmap();
        scene->setDirty();
    }
}

bool OptimizationContext::loadMesh(std::string &path){
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
        // @todo synchronization
        m_meshData = importer->mesh(0);
        if(!m_meshData) return false;

        std::tie(m_V,m_F) = toEigen(*m_meshData);
        m_cost->emplace_back<InterfaceEnergy>("modical mortola", m_weightMM, m_V, m_F, m_epsilon);
        m_cost->emplace_back<PotentialEnergy>("modical mortola", m_weightMM, m_V, m_F, m_epsilon);
        m_cost->emplace_back<AreaRegularizer>("area regularizer", m_weightArea, m_V, m_F);
        m_cost->emplace_back<ConnectednessConstraint>("connectedness constraint", m_weightConnectedness, m_V, m_F, m_epsilon, m_a, m_b);
        m_cost->emplace_back<ConnectednessConstraint>("connectedness constraint", m_weightConnectedness, m_V, m_F, m_epsilon, -m_b, -m_a);
        return true;
    }
    else{
        Debug{} << "Could not load mesh";
        return false;
    }
}

template <class To, class From>
auto bit_cast(const From &src) noexcept
{
    To dst;
    std::memcpy(&dst, &src, sizeof(To));
    return dst;
}

bool OptimizationContext::saveMesh(std::string &path) {
    //@todo generate path if not existent


    auto* mesh = new aiMesh;
    if(m_meshData->isIndexed()){
        CORRADE_INTERNAL_ASSERT(m_meshData->primitive() == MeshPrimitive::Triangles);
        auto numTriangles = m_meshData->indexCount() / 3;
        mesh->mNumFaces = numTriangles;
        auto faces = new aiFace[numTriangles];
        for (int i = 0; i < numTriangles; ++i) {
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

    std::transform(vertexView.begin(), vertexView.end(), vertices, bit_cast<aiVector3D, Vector3>);
    std::transform(normalView.begin(), normalView.end(), normals, bit_cast<aiVector3D, Vector3>);
    std::transform(colorView.begin(), colorView.end(), colors, bit_cast<aiColor4D, Color4>);

    mesh->mColors[0] = colors;
    mesh->mVertices = vertices;
    mesh->mNormals = normals;

    aiScene scene;
    auto meshes = new aiMesh*;
    meshes[0] = mesh;
    scene.mMeshes = meshes;

    Assimp::Exporter exporter;
    auto desc = exporter.GetExportFormatDescription(0);
    exporter.Export(&scene, desc->id, path.c_str(), 0);
}

void OptimizationContext::showMenu(ImGuiIntegration::Context& context){
    ImGui::SetNextWindowPos({500.0f, 50.0f}, ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowBgAlpha(0.5f);
    ImGui::Begin("Options", nullptr);
    ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.6f);

    /* General information */
    ImGui::Text("Hide/show menu: H");
    ImGui::Text("Rendering: %3.2f FPS (1 thread)", Double(ImGui::GetIO().Framerate));
    ImGui::Spacing();

    /* Rendering parameters */
    if(ImGui::TreeNodeEx("Particle Rendering", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::PushID("Particle Rendering");
        {
            constexpr const char* items[] = {"Uniform", "Ramp by ID", "Random"};
            static Int colorMode = 1;
            if(ImGui::Combo("Color Mode", &colorMode, items, 3));

            if(colorMode == 0) { /* Uniform color */
                static Color3 color = Color3::blue();
                if(ImGui::ColorEdit3("Diffuse Color", color.data()));
            }
        }
        static Vector3 lightDir = Vector3::xAxis();
        if(ImGui::InputFloat3("Light Direction", lightDir.data()));
        ImGui::PopID();
        ImGui::TreePop();
    }
    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    /* Simulation parameters */
    if(ImGui::TreeNodeEx("Simulation", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::PushID("Simulation");
        Float a,b,c;
        bool d;
        ImGui::InputFloat("Stiffness", &a);
        ImGui::SliderFloat("Viscosity", &b, 0.0f, 1.0f);
        ImGui::SliderFloat("Restitution", &c, 0.0f, 1.0f);
        ImGui::Checkbox("Dynamic Boundary", &d);
        ImGui::PopID();
        ImGui::TreePop();
    }
    ImGui::Spacing();
    ImGui::Separator();

    /* Reset */
    ImGui::Spacing();
    bool p = false;
    if(ImGui::Button(p ? "Play Sim" : "Pause Sim"))
        p ^= true;
    ImGui::SameLine();
    if(ImGui::Button("Reset Sim")) {
        p = false;
    }
    ImGui::SameLine();
    if(ImGui::Button("Reset Camera"));
    ImGui::PopItemWidth();
    ImGui::End();
}

void OptimizationContext::work()
{
    while(m_exit){
        std::unique_lock ul(m_mutex);
        if(m_run){
           ul.unlock();
           Solve::Summary summary;
           Solve(m_options, m_problem, m_U, summary);
           continue;
        }
        else
           m_cv.wait(ul, [this]{ return m_run; });
    }
}

OptimizationContext::~OptimizationContext(){
    m_exit = true;

    {
        std::lock_guard l(m_mutex);
        m_run = true;
    }

    m_cv.notify_all();

    Debug{} << "Waiting for optimization thread to join...";
    m_optThread.join();
    Debug{} << "Done!";
}
