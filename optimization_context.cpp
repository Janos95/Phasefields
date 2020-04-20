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

OptimizationContext::OptimizationContext(PhasefieldData& data):
    pd(data),
    options{
        .max_num_iterations = 10000,
        .minimizer_progress_to_stdout = true,
        .update_state_every_iteration = true,
        .callbacks = {[this]{}}
             }
{
}

solver::Status OptimizationContext::optimizationCallback(Solver::IterationSummary const&){
    if(!optimize)
        return Solver::Status::Abort;

    {
        std::lock_guard l(pd.mutex);
        auto textureCoordsView = m_pd.meshData.mutableAttribute<Vector2>(Trade::MeshAttribute::TextureCoordinates);
        CORRADE_INTERNAL_ASSERT(textureCoordsView.size() == m_U.size());
        for (int i = 0; i < textureCoordsView.size(); ++i)
            textureCoordsView[i].x() = static_cast<float>(.5 * (m_U[i] + 1.));

        if(m_showShortestPaths){

        }

        pd.status = PhasefieldData::Status::PhasefieldUpdated;
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
    taskGroup.wait();

    optimize = true;
    taskGroup.run([this]{
        Solver::Summary summary;
        //ceres::GradientProblemSolver::Summary summary;
        //ceres::GradientProblem problem(m_cost);
        //ceres::Solve(options, problem, m_U.data(), &summary);
        Solver(options, problem, phasefield.data(), summary);
        Cr::Utility::copy(phasefield, pd.phasefield);
        std::copy(phas.begin(), m_U.end(), m_pd.phasefield.begin());
        optimize = false;
        Debug{} << summary.briefReport().c_str();
            });
}

void OptimizationContext::stopOptimization() {
    optimize = false;
}

void drawConnectednessConstraintOptions(ConnectednessConstraint<Float>& f){

    ImGui::DragFloatRange2("Positive Interval", &f.a, &f.b, .01f, .0f, 1.f, "Min: %.2f", "Max: %.2f");
}

void drawDirichletEnergyOptions(DirichletEnergy& f){

}

void drawDoubleWellPotentialOptions(DoubleWellPotential& f){

}

void drawAreaRegularizationOptions(AreaRegularizer& f){

}

void OptimizationContext::drawImGui(Viewer&){
    if (ImGui::TreeNode("Optimization Options"))
    {
        static std::uint32_t iterations = 100;
        constexpr static std::uint32_t step = 1;
        ImGui::InputScalar("iterations", ImGuiDataType_S32, &m_options.max_num_iterations, &step, nullptr, "%u");

        constexpr float minEps = 0.f, maxEps = 1.;
        ImGui::DragScalar("epsilon", ImGuiDataType_Double, &m_epsilon, .01f, &minEps, &maxEps, "%f", 2);

        ImGui::Text("Preimage Intervals");
        ImGui::DragFloatRange2("Negative Interval", &m_negA, &m_negB, .01f, -1.f, .0f, "Min: %.2f", "Max: %.2f");

        ImGui::Text("Functional Weights");
        for(auto& [name, w] : m_weights){
            constexpr static double lower = 0., upper = 1.;
            ImGui::SliderScalar(name.c_str(), ImGuiDataType_Double, &w, &lower, &upper, "%.5f", 1.0f);
        }

        for(auto& f : pd.functionals){
            switch(f->type){
                case FunctionalType::DoubleWellPotential :
                case FunctionalType::DirichletEnergy :
                case FunctionalType::Area : draw
                case FunctionalType::Connectedness : drawConnectednessConstraintOptions(dynamic_cast<ConnectednessConstraint<Float>&>(*f));
            }
        }

        if (ImGui::Button("Optimize"))
            startOptimization();

        if (auto stop = ImGui::Button("Stop"))
            stopOptimization();

        ImGui::TreePop();
    }
}

void OptimizationContext::tickEvent(Viewer& viewer){

}
