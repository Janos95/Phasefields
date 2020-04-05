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
#include <Magnum/EigenIntegration/Integration.h>
#include <Magnum/ImGuiIntegration/Context.hpp>

#include <folly/futures/Future.h>
#include <folly/executors/ThreadedExecutor.h>
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

void OptimizationContext::startOptimization() {

    folly::Future<folly::Unit> m_optimizationFuture;
    folly::ThreadedExecutor m_executor;

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


/**
 *
 * @todo
 *  - add color edit for phasefield vis
 *  - add options to change epsilon
 */
void OptimizationContext::drawEvent(){
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


    ImGui::Checkbox("Enable Brush", &m_brushEnabled);

    ImGui::PopItemWidth();
    ImGui::End();
}

OptimizationContext::~OptimizationContext(){
    if(m_optimizationFuture.valid() && !m_optimizationFuture.isReady()){
        m_exit = true;
        m_optimizationFuture.wait();
    }
}
