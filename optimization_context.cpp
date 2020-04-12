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

#include <Magnum/EigenIntegration/Integration.h>
#include <Magnum/ImGuiIntegration/Context.hpp>


using namespace Magnum;
using namespace Corrade;


void toEigen(PhasefieldData& pd, Eigen::MatrixXi& F, Eigen::MatrixXd& V, Eigen::VectorXd& U){
    auto& vertices = pd.V;
    auto& indices = pd.F;
    auto& phasefield = pd.phasefield;

    using MatrixXU = Eigen::Matrix<UnsignedInt, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    F = Eigen::Map<const MatrixXU>(indices.data(), indices.size() / 3, 3).cast<int>();
    V = Eigen::MatrixXd(vertices.size(), 3);
    for (std::size_t i = 0; i < pd.V.size(); ++i)
        V.row(i) = EigenIntegration::cast<Eigen::Vector3f>(vertices[i]).cast<double>();

    U = Eigen::Map<Eigen::VectorXf>(phasefield.data(), phasefield.size()).cast<double>();
}


OptimizationContext::OptimizationContext(PhasefieldData& data):
    m_pd(data),
    m_optimizationCallback([this](auto const& summary){ return optimizationCallback(summary); }),
    m_options{
        .max_num_iterations = 100,
        .minimizer_progress_to_stdout = true,
        .update_state_every_iteration = true,
        .callbacks = {&m_optimizationCallback}
             },
    m_cost(new SumProblem())
{
}

ceres::CallbackReturnType OptimizationContext::optimizationCallback(ceres::IterationSummary const&){
    if(!m_continue)
        return ceres::SOLVER_TERMINATE_SUCCESSFULLY;

    {
        std::lock_guard l(m_pd.mutex);
        auto textureCoordsView = m_pd.meshData.mutableAttribute<Vector2>(Trade::MeshAttribute::TextureCoordinates);
        CORRADE_INTERNAL_ASSERT(textureCoordsView.size() == m_U.size());
        for (int i = 0; i < textureCoordsView.size(); ++i)
            textureCoordsView[i].x() = static_cast<float>(.5 * (m_U[i] + 1.));
        m_pd.status = PhasefieldData::Status::PhasefieldUpdated;
    }

    return ceres::SOLVER_CONTINUE;
}

#define EMPLACE_COST(cost, func, eps, ...) if(eps > std::numeric_limits<double>::epsilon()) cost->emplace_back<func>(__VA_ARGS__);

/*
 * @todo only need to construct the cost functional when we get a new mesh/subdivision event.
 * But then we would need to handle the resetting of parameters which at the moment is a bit combersome...
 */
void OptimizationContext::startOptimization() {

    stopOptimization();

    toEigen(m_pd, m_F, m_V, m_U);
    m_cost->clear();

    if(m_weights[0].second > std::numeric_limits<double>::epsilon())
        m_cost->emplace_back<InterfaceEnergy>("Modica Mortola", 1., m_V, m_F, m_epsilon);
    if(m_weights[0].second > std::numeric_limits<double>::epsilon())
        m_cost->emplace_back<PotentialEnergy>("Modica Mortola", 1., m_V, m_F, m_epsilon);
    if(m_weights[1].second > std::numeric_limits<double>::epsilon())
        m_cost->emplace_back<AreaRegularizer>("Area Regularization", 1., m_V, m_F);
    if(m_weights[2].second > std::numeric_limits<double>::epsilon())
        m_cost->emplace_back<ConnectednessConstraint<double>>("Connectedness Constraint Positive", 1., m_V, m_F, m_epsilon, m_posA, m_posB);
    if(m_weights[3].second > std::numeric_limits<double>::epsilon())
        m_cost->emplace_back<ConnectednessConstraint<double>>("Connectedness Constraint Negative", 1., m_V, m_F, m_epsilon, -m_negA, -m_negB);

    if(!m_problem)
        m_problem.emplace(m_cost);

    m_continue = true;
    m_thread = std::thread([&, options = m_options]{
        ceres::GradientProblemSolver::Summary summary;
        ceres::Solve(options, *m_problem, m_U.data(), &summary);
        std::copy(m_U.begin(), m_U.end(), m_pd.phasefield.begin());
            });
}

void OptimizationContext::stopOptimization() {
    m_continue = false;
    if(m_thread.joinable())
        m_thread.join();
}


void OptimizationContext::drawImGui(){
    if (ImGui::TreeNode("Optimization Options"))
    {
        static std::uint32_t iterations = 100;
        constexpr static std::uint32_t step = 1;
        ImGui::InputScalar("iterations", ImGuiDataType_S32, &m_options.max_num_iterations, &step, nullptr, "%u");

        constexpr float minEps = 0.f, maxEps = 1.;
        ImGui::DragScalar("epsilon", ImGuiDataType_Double, &m_epsilon, .01f, &minEps, &maxEps, "%f", 2);

        ImGui::Text("Preimage Intervals");
        ImGui::DragFloatRange2("Positive Interval", &m_posA, &m_posB, .01f, .0f, 1.f, "Min: %.2f", "Max: %.2f");
        ImGui::DragFloatRange2("Negative Interval", &m_negA, &m_negB, .01f, -1.f, .0f, "Min: %.2f", "Max: %.2f");

        ImGui::Text("Functional Weights");
        for(auto& [name, w] : m_weights){
            constexpr static double lower = 0., upper = 1.;
            ImGui::SliderScalar(name.c_str(), ImGuiDataType_Double, &w, &lower, &upper, "%.5f", 1.0f);
        }

        if (ImGui::Button("Optimize"))
            startOptimization();

        if (auto stop = ImGui::Button("Stop"))
            stopOptimization();

        ImGui::TreePop();
    }
}

OptimizationContext::~OptimizationContext(){
    stopOptimization();
}
