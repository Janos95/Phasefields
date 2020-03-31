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

#include <random>

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
    m_optimizationCallback([this](auto const& summary){ return optimizationCallback(summary); }),
    m_options{
        .max_num_iterations = 100,
        .minimizer_progress_to_stdout = true,
        .update_state_every_iteration = true,
        .callbacks = {&m_optimizationCallback}},
    m_cost(new SumProblem()),
    m_inputPath(100),
    m_outputPath(100)
{
    std::string defaultPath = "/home/janos/data/spot.ply";
    std::copy(defaultPath.begin(), defaultPath.end(), m_inputPath.begin());
}

ceres::CallbackReturnType OptimizationContext::optimizationCallback(ceres::IterationSummary const&){
    bool exit = true;
    if(m_exit.compare_exchange_strong(exit, false)){
        return ceres::SOLVER_TERMINATE_SUCCESSFULLY;
    }

    {
        std::lock_guard l(m_mutex);
        auto colorView = m_meshData->mutableAttribute<Color4>(Trade::MeshAttribute::Color);
        std::transform(m_U.begin(), m_U.end(), colorView.begin(), JetColorMap{});
    }

    m_reupload = true;
    return ceres::SOLVER_CONTINUE;
}

void OptimizationContext::updateScene(Scene*& scene){
    bool newMeshData = true;
    if(m_newMeshData.compare_exchange_strong(newMeshData, false)) {
        Debug{} << "Uploading new mesh";
        if(scene) scene->reset();
        else scene = new Scene(); //@todo scene leaks, but I am not sure who should own it...
        {
            std::lock_guard l(m_mutex);
            scene->addObject("mesh", upload(*m_meshData));
        }
        scene->addNode("mesh", "mesh", ShaderType::Phong);
        scene->setDirty();
    }
    bool reupload = true;
    if(m_reupload.compare_exchange_strong(reupload, false)){
        auto* obj = scene->getObject("mesh");
        CORRADE_INTERNAL_ASSERT(obj);
        GL::Buffer& vertices = obj->vertices;
        auto data = vertices.map(0,
                                 vertices.size(),
                                 GL::Buffer::MapFlag::Write);

        CORRADE_CONSTEXPR_ASSERT(data, "could not map vertex data");
        {
            std::lock_guard l(m_mutex);
            Utility::copy(m_meshData->vertexData(), data);
        }
        vertices.unmap();
        scene->setDirty();
    }
}

void OptimizationContext::initializePhasefield(InitializationFlag flag)
{
    std::default_random_engine engine(0);
    std::normal_distribution distr(.0, .1);

    if(m_optimizationFuture.valid() && !m_optimizationFuture.isReady()){
        m_exit = true;
        m_optimizationFuture.wait();
    }

    for(auto& u: m_U) u = distr(engine);
    auto colorView = m_meshData->mutableAttribute<Color4>(Trade::MeshAttribute::Color);
    std::transform(m_U.begin(), m_U.end(), colorView.begin(), JetColorMap{});
    m_reupload = true;
}

bool OptimizationContext::loadMesh(std::string const& path){
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


void OptimizationContext::startOptimization() {
    stopOptimization();

    if(!m_problem)
        return;

    fmt::print("Optimizing the following functionals for "
               "a maximum of {} iterations: \n", m_options.max_num_iterations);
    for(auto const& [name, w] : m_checkBoxes){
        m_cost->setWeight(name, w);
        if(w){
            fmt::print("{}\n",name);
        }
    }

    m_cost->visit("Connectedness Constraint Positive",[=](auto& p ){ auto ptr = p.get(); dynamic_cast<ConnectednessConstraint*>(ptr)->setPreimageInterval(m_posA, m_posB); });
    m_cost->visit("Connectedness Constraint Negative",[=](auto& p ){ auto ptr = p.get(); dynamic_cast<ConnectednessConstraint*>(ptr)->setPreimageInterval(m_negA, m_negB); });

    m_optimizationFuture = folly::makeSemiFuture().via(&m_executor).thenValue(
            [this](auto&&){
                ceres::GradientProblemSolver::Options options;
                {
                    std::lock_guard l(m_mutex);
                    options = m_options;
                }

                ceres::GradientProblemSolver::Summary summary;
                ceres::Solve(options, *m_problem, m_U.data(), &summary);
            });
}

void OptimizationContext::stopOptimization() {
    if(m_optimizationFuture.valid() && !m_optimizationFuture.isReady()){
        m_exit = true;
        m_optimizationFuture.wait();
    }
}

template <class To, class From>
To bit_cast(const From &src) noexcept
{
    To dst;
    std::memcpy(&dst, &src, sizeof(To));
    return dst;
}

bool OptimizationContext::saveMesh(std::string const& path) {
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

/**
 *
 * @todo
 *  - add color edit for phasefield vis
 *  - add save option
 *  - add options to change parameters
 *  - add options to enable/disable functionals
 */
void OptimizationContext::showMenu(ImGuiIntegration::Context& context){
    ImGui::SetNextWindowPos({500.0f, 50.0f}, ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowBgAlpha(0.5f);
    ImGui::Begin("Options", nullptr);
    ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.6f);


    ImGui::InputText("mesh path", m_inputPath.data(), m_inputPath.size());
    if(ImGui::Button("load mesh")){
        loadMesh(std::string{m_inputPath.begin(),std::find(m_inputPath.begin(), m_inputPath.end(), '\0')});
    }


    if(ImGui::Button("reset phasefield")){
        initializePhasefield(InitializationFlag::RandomNormal);
    }

    if (ImGui::Button("optimize")){
        startOptimization();
    }

    static std::uint32_t iterations = 100;
    constexpr static std::uint32_t step = 1;
    std::uint32_t old = iterations;
    ImGui::InputScalar("iterations", ImGuiDataType_U32, &iterations, &step, nullptr, "%u");
    if(old != iterations){
        std::lock_guard l(m_mutex);
        m_options.max_num_iterations = static_cast<int>(iterations);
    }

    for(auto& [name, w] : m_checkBoxes){
        constexpr static double lower = 0., upper = 1.;
        ImGui::SliderScalar(name.c_str(), ImGuiDataType_Double, &w, &lower, &upper, "%.5f", 1.0f);
    }

    if (ImGui::TreeNode("Preimage Intervals"))
    {
        ImGui::DragFloatRange2("Positive Interval", &m_posA, &m_posB, .01f, .0f, 1.f, "Min: %.2f", "Max: %.2f");
        ImGui::DragFloatRange2("Negative Interval", &m_negA, &m_negB, .01f, -1.f, .0f, "Min: %.2f", "Max: %.2f");
        ImGui::TreePop();
    }

    ImGui::PopItemWidth();
    ImGui::End();
}

OptimizationContext::~OptimizationContext(){
    if(m_optimizationFuture.valid() && !m_optimizationFuture.isReady()){
        m_exit = true;
        m_optimizationFuture.wait();
    }
}
