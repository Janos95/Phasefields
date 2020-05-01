//
// Created by janos on 26.03.20.
//

#include "optimization_context.hpp"
#include "modica_mortola.hpp"
#include "connectedness_constraint.hpp"
#include "upload.hpp"
#include "loss_functions.hpp"
#include "custom_widgets.hpp"

#include <Corrade/Containers/Pointer.h>
#include <Corrade/Utility/Algorithms.h>

#include <Magnum/MeshTools/Compile.h>
#include <Magnum/ImGuiIntegration/Context.hpp>
#include <imgui_internal.h>

#include <iostream>

using namespace Magnum;
using namespace Corrade;

namespace {

struct ComboElement {
    std::string name;
    LossType type;
};


struct ComboElementFunctional {
    std::string name;
    FunctionalType type;
};

struct CompareByName {
    bool operator()(ComboElement const &e1, ComboElement const &e2) {
        return e1.name < e2.name;
    }
};

auto makeComboMapFunctionals() {
    Containers::Array<ComboElementFunctional> map{Containers::InPlaceInit,
        {
            {"Dirichlet Energy", FunctionalType::DirichletEnergy},
            {"Double Well Potential", FunctionalType::DoubleWellPotential},
            {"Area Quadratic", FunctionalType::Area1},
            {"Area Smooth Step", FunctionalType::Area2},
            {"Connectedness Constraint", FunctionalType::Connectedness},
        }};
    return map;
}




auto makeComboMap() {
    Containers::Array<ComboElement> map{Containers::InPlaceInit,
        {
            {"Composed Loss", LossType::Composed},
            {"Scaled Loss", LossType::Scaled},
            {"Tolerant Loss", LossType::Tolerant},
            {"Trivial Loss", LossType::Trivial},
            {"Huber Loss", LossType::Huber},
            {"Soft One Loss", LossType::SoftOne},
            {"Cauchy Loss", LossType::Cauchy},
            {"Arctan Loss", LossType::Arctan},
            {"Tukey Loss", LossType::Tukey},
            {"Quadratic Loss", LossType::Quadratic}
        }};
    return map;
}

}


Path::Path(Object3D* parent, Trade::MeshData const& md) :
        Object3D(parent),
        mesh(MeshTools::compile(md))
{
}


template<class Scalar>
struct UpdateVis {
    OptimizationContext& context;
    ConnectednessConstraint<Scalar>& connectednessConstraint;

    bool drawShortestPaths = false;
    bool drawConnectedComponents = false;
    bool drawGeodesicDistance = false;

    Cr::Containers::Array<int> components;
    Cr::Containers::Array<Double> geodesicDistance;
    int numComponents;

    solver::Status operator()(solver::IterationSummary const&) {
        std::lock_guard l(context.mutex);
        if(!context.optimize) return solver::Status::ABORT;
        if(drawShortestPaths)
            context.paths = std::move(connectednessConstraint.lineStrips);
        if(drawConnectedComponents){
            Utility::copy(connectednessConstraint.components, components);
            numComponents = connectednessConstraint.numComponents;
        } else if(drawGeodesicDistance){
            Utility::copy(connectednessConstraint.components, geodesicDistance);
        }
        context.update = true;
        return solver::Status::CONTINUE;
    }
};

using uf = unique_function<solver::Status(solver::IterationSummary const&)>;

OptimizationContext::OptimizationContext(PhasefieldData& data):
    pd(data),
    options{
        .max_num_iterations = 10000,
        .minimizer_progress_to_stdout = true,
        .update_state_every_iteration = true,
    },
    dirichletScaling(1.), doubleWellScaling(1.), connectednessScaling(1.)
{
    Containers::arrayAppend(options.callbacks, Containers::InPlaceInit, [this](auto const& s){ return updatePhasefield(s); });
}

void OptimizationContext::startOptimization() {

    stopOptimization();

    Containers::arrayResize(phasefield, pd.phasefield.size());
    Containers::arrayResize(parameters, pd.phasefield.size());
    Utility::copy(pd.phasefield, parameters);

    optimize = true;
    thread = std::thread([this]{
        solver::Summary summary;
        solve(options, pd.problem, parameters.data(), &summary);
        //Debug{} << summary.briefReport().c_str();
    });
}

void OptimizationContext::stopOptimization() {
    if(optimize.exchange(false))
        thread.join();
}

void drawLoss(Containers::Pointer<LossFunction>& loss, int& nodeCount){
    if(loss->type == LossType::Unknown) return;
    ImGui::PushID(nodeCount);
    CORRADE_ASSERT(loss, "Loss is null",);
    //Chose a loss
    static auto map = makeComboMap();
    ComboElement* cur = std::find_if(map.begin(), map.end(), [&](auto& el){ return el.type == loss->type; });
    bool newLoss = false;
    if (ImGui::BeginCombo("##combo", cur->name.c_str(), ImGuiComboFlags_NoArrowButton)){
        for (auto it = map.begin(); it < map.end(); ++it){
            bool isSelected = (cur == it);
            if (ImGui::Selectable(it->name.c_str(), isSelected)){
                cur = it;
                newLoss = true;
            }
            if (isSelected)
                ImGui::SetItemDefaultFocus();
        }
        ImGui::EndCombo();
    }

    constexpr Double min = 0.f, max = 1.;
    switch (cur->type) {
        case LossType::Tolerant : {
            if(newLoss) loss = Containers::pointer<TolerantLoss>(1., 1.);
            auto& l = dynamic_cast<TolerantLoss&>(*loss);
            auto updateC = ImGui::DragScalar("a", ImGuiDataType_Double, &l.a, .01f, &min, &max, "%f", 2);
            updateC |= ImGui::DragScalar("b", ImGuiDataType_Double, &l.b, .01f, &min, &max, "%f", 2);
            if(updateC) l.c =  l.b * log(1.0 + exp(-l.a / l.b));
            break;
        }
        case LossType::Trivial : {
            if(newLoss) loss = Containers::pointer<TrivialLoss>();
            break;
        }
        case LossType::Quadratic : {
            if(newLoss) loss = Containers::pointer<QuadrationLoss>();
            break;
        }
        case LossType::Composed : {
            if(newLoss) loss = Containers::pointer<ComposedLoss>(Containers::pointer<TrivialLoss>(), Containers::pointer<TrivialLoss>());
            auto& l = dynamic_cast<ComposedLoss&>(*loss);
            ImGui::Indent();
            drawLoss(l.f, ++nodeCount);
            drawLoss(l.g, ++nodeCount);
            ImGui::Unindent();
            break;
        }
        case LossType::Scaled : {
            if(newLoss) loss = Containers::pointer<ScaledLoss>(Containers::pointer<TrivialLoss>(), 1.);
            auto& l = dynamic_cast<ScaledLoss&>(*loss);
            ImGui::DragScalar("Weight", ImGuiDataType_Double, &l.a, .01f, &min, &max, "%f", 2);
            ImGui::Indent();
            drawLoss(l.f, ++nodeCount);
            ImGui::Unindent();
            break;
        }
        case LossType::Cauchy : {
            if(newLoss) loss = Containers::pointer<CauchyLoss>(1.);
            auto& l = dynamic_cast<CauchyLoss&>(*loss);
            ImGui::SliderScalar("b", ImGuiDataType_Float, &l.b, &min, &max, "%f");
            break;
        }
        case LossType::Arctan : {
            if(newLoss) loss = Containers::pointer<ArctanLoss>(1.);
            auto& l = dynamic_cast<ArctanLoss&>(*loss);
            break;
        }
        case LossType::Tukey : {
            if(newLoss) loss = Containers::pointer<TukeyLoss>(1.);
            auto& l = dynamic_cast<TukeyLoss&>(*loss);
            break;
        }
        case LossType::SoftOne : {
            if(newLoss) loss = Containers::pointer<SoftLOneLoss>(1.);
            auto& l = dynamic_cast<SoftLOneLoss&>(*loss);
            break;
        }
        case LossType::Huber : {
            if(newLoss) loss = Containers::pointer<HuberLoss>(1.);
            auto& l = dynamic_cast<HuberLoss&>(*loss);
            break;
        }
        case LossType::Unknown : assert(false);
    }
    ImGui::PopID();
}



void OptimizationContext::drawConnectednessConstraintOptions(ConnectednessConstraint<Double>& f){
    ImGui::Text("Connectedness Constraint");
    dragDoubleRange2("Positive Interval", &f.a, &f.b, 0.01f, .0f, 1.f, "Min: %.2f", "Max: %.2f",1.f);
    ImGui::Text("Visualize Shortest Paths");
    ImGui::SameLine();
    toggleButton("##button", &drawShortestPaths);
}

void drawDirichletEnergyOptions(DirichletEnergy& f){
    ImGui::Text("Dirichlet Energy");
}

void drawDoubleWellPotentialOptions(DoubleWellPotential& f){
    ImGui::Text("Double Well Potential");
}

void drawAreaRegularization1Options(AreaRegularizer1& f){
    ImGui::Text("Area Regularization");
    constexpr double min = 0.f, max = 1.;
    ImGui::SliderScalar("Area Ratio", ImGuiDataType_Double, &f.areaRatio, &min, &max, "%.5f", 1.0f);
}


void drawAreaRegularization2Options(AreaRegularizer2& f){
    ImGui::Text("Area Regularization");
    constexpr double min = 0.f, max = 1.;
    ImGui::SliderScalar("Area Ratio", ImGuiDataType_Double, &f.areaRatio, &min, &max, "%.5f", 1.0f);
}

void OptimizationContext::appendFunctional(Containers::Array<Containers::Pointer<Functional>>& functionals, FunctionalType type){
    auto triangles = Containers::arrayCast<Vector3ui>(pd.F);
    switch (type) {
        case FunctionalType::Area1 : {
            auto p = Containers::pointer<AreaRegularizer1>(pd.V, triangles);
            p->type = FunctionalType::Area1;
            Containers::arrayAppend(functionals, Containers::InPlaceInit, std::move(p));
            break;
        }
        case FunctionalType::Area2 : {
            auto p = Containers::pointer<AreaRegularizer2>(pd.V, triangles);
            p->type = FunctionalType::Area2;
            Containers::arrayAppend(functionals, Containers::InPlaceInit, std::move(p));
            break;
        }
        case FunctionalType::DirichletEnergy : {
            auto p = Containers::pointer<DirichletEnergy>(pd.V, triangles);
            p->scaling = dirichletScaling;
            Containers::arrayAppend(functionals, Containers::InPlaceInit, std::move(p));
            break;
        }
        case FunctionalType::DoubleWellPotential : {
            auto p = Containers::pointer<DoubleWellPotential>(pd.V, triangles);
            p->type = FunctionalType::DoubleWellPotential;
            p->scaling = doubleWellScaling;
            Containers::arrayAppend(functionals, Containers::InPlaceInit, std::move(p));
            break;
        }
        case FunctionalType::Connectedness :
            auto p = Containers::pointer<ConnectednessConstraint<Double>>(pd.V, triangles, .1f, .95f);
            p->scaling = connectednessScaling;
            p->drawLineStrips = &drawShortestPaths;
            Containers::arrayAppend(options.callbacks, Containers::InPlaceInit, UpdateVis<Double>{*this, *p});
            Containers::arrayAppend(functionals, Containers::InPlaceInit, std::move(p));
            break;
    }
}

void OptimizationContext::drawImGui(Viewer&){
    if (ImGui::TreeNode("Optimization Options"))
    {
        static auto map = makeComboMapFunctionals();
        static auto cur = map.end();
        if(ImGui::BeginCombo("##combo", cur != map.end() ? cur->name.c_str() : nullptr)){
            for (auto it = map.begin(); it < map.end(); ++it){
                bool isSelected = (cur == it);
                if (ImGui::Selectable(it->name.c_str(), isSelected)){
                    cur = it;
                }
                if (isSelected)
                    ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }

        if(ImGui::Button("Add Functional As Cost") && cur != map.end())
            appendFunctional(pd.problem.functionals, cur->type);

        ImGui::SameLine();

        if(ImGui::Button("Add Functional As Constraint") && cur != map.end())
            appendFunctional(pd.problem.constraints, cur->type);


        int nodeCount = 0;
        int toRemove = -1;
        auto& fs = pd.problem.functionals;
        for(int i = 0; i < fs.size(); ++i){
            auto& f = *fs[i];
            ImGui::PushID(++nodeCount);
            ImGui::Separator();
            if(ImGui::Button("Remove Functional")){
                toRemove = i;
            }
            drawLoss(f.loss, nodeCount);
            switch(f.type){
                case FunctionalType::DoubleWellPotential :
                    drawDoubleWellPotentialOptions(dynamic_cast<DoubleWellPotential&>(f));
                    break;
                case FunctionalType::DirichletEnergy :
                    drawDirichletEnergyOptions(dynamic_cast<DirichletEnergy&>(f));
                    break;
                case FunctionalType::Area1 :
                    drawAreaRegularization1Options(dynamic_cast<AreaRegularizer1&>(f));
                    break;
                case FunctionalType::Area2 :
                    drawAreaRegularization2Options(dynamic_cast<AreaRegularizer2&>(f));
                    break;
                case FunctionalType::Connectedness :
                    drawConnectednessConstraintOptions(dynamic_cast<ConnectednessConstraint<Double>&>(f));
                    break;
            }
            ImGui::PopID();
        }

        if(toRemove >= 0){
            std::swap(fs[toRemove], fs.back());
            Containers::arrayResize(fs, fs.size() - 1);
        }

        static Double epsilon = 0.1;
        if(dirichletScaling.refCount() > 1 || doubleWellScaling.refCount() > 1 || connectednessScaling.refCount() > 1){
            constexpr Mg::Double minEps = 0.f, maxEps = 1.;
            ImGui::DragScalar("epsilon", ImGuiDataType_Double, &epsilon, .01f, &minEps, &maxEps, "%f", 2);
            *dirichletScaling = epsilon/2.;
            *doubleWellScaling = 1./epsilon;
            *connectednessScaling = 1./(epsilon * epsilon);
        }


#ifndef NDEBUG
        if(ImGui::Button("debug button")) {
            std::size_t n = pd.problem.numParameters();
            double r;
            Containers::Array<Double> numericGrad(Containers::NoInit, n);
            Containers::Array<Double> grad(Containers::NoInit, n);
            Containers::Array<Double> params(Containers::NoInit, n);
            Utility::copy(pd.phasefield, params);
            for (auto& f : pd.problem.functionals) {
                std::cout << " ================================ " << std::endl;
                f->evaluate(pd.phasefield.data(), &r, grad.data());
                for (int j = 0; j < n; ++j) {
                    double r1, r2;
                    double h = 10e-3;
                    params[j] = pd.phasefield[j] + h;
                    f->evaluate(params.data(), &r1, nullptr);
                    params[j] = pd.phasefield[j] - h;
                    f->evaluate(params.data(), &r2, nullptr);
                    numericGrad[j] = (r1 - r2) / (2 * h);
                    params[j] = pd.phasefield[j];
                }
                Debug{} << grad;
                Debug{} << numericGrad;
            }

            std::cout << " ================================ " << std::endl;
            std::cout << " whole problem" << std::endl;
            pd.problem.evaluate(pd.phasefield.data(), &r, grad.data());
            for (int j = 0; j < n; ++j) {
                double r1, r2;
                double h = 10e-3;
                params[j] = pd.phasefield[j] + h;
                pd.problem.evaluate(params.data(), &r1, nullptr);
                params[j] = pd.phasefield[j] - h;
                pd.problem.evaluate(params.data(), &r2, nullptr);
                numericGrad[j] = (r1 - r2) / (2 * h);
                params[j] = pd.phasefield[j];
            }
            Debug{} << grad;
            Debug{} << numericGrad;

        }
#endif
        static std::uint32_t iterations = 100;
        constexpr static std::uint32_t step = 1;
        ImGui::InputScalar("iterations", ImGuiDataType_S32, &options.max_num_iterations, &step, nullptr, "%u");

        if (ImGui::Button("Optimize") && !pd.problem.functionals.empty())
            startOptimization();

        ImGui::SameLine();

        if (auto stop = ImGui::Button("Stop"))
            stopOptimization();

        ImGui::TreePop();
    }
}

void tickEvent(Viewer& viewer){
    if(update.exchange(false)){

        if(drawShortestPaths){
            auto& manipulator = *viewer.pathManipulator;
            manipulator.children().clear();
            std::lock_guard l(mutex);
            for(auto const& md : paths){
                auto p = new Path(&manipulator, md);
                auto drawable = new ColorMapFlatDrawable(*p, p->mesh, pathShader, &viewer.drawableGroup);
                drawable->color = pathColor;
            }
        }
        auto texture = pd.meshData.mutableAttribute<Vector2>(Trade::MeshAttribute::TextureCoordinates);
        {
            std::lock_guard l(mutex);
            if(drawConnectedComponents){
                auto tf = 1.f/static_cast<float>(numComponents);
                for (int i = 0; i < phasefield.size(); ++i)
                    texture[i].x() = tf * components[i];
            } else {
                for (int i = 0; i < phasefield.size(); ++i)
                    texture[i].x() = .5f * (pd.phasefield[i] + 1.f);
            }
        }

        reuploadVertices(pd, pd.meshData);
        pd.status = PhasefieldData::Status::OptimizationUpdated;
    } else if(pd.status == PhasefieldData::Status::OptimizationUpdated) {
        pd.status == PhasefieldData::Status::NothingChanged;
    }
}
