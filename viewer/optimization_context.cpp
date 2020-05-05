//
// Created by janos on 26.03.20.
//

#include "optimization_context.hpp"
#include "modica_mortola.hpp"
#include "upload.hpp"
#include "loss_functions.hpp"
#include "viewer.hpp"

#include <Corrade/Containers/Pointer.h>

#include <Magnum/ImGuiIntegration/Context.hpp>


using namespace Magnum;
using namespace Corrade;

namespace {

struct ComboElement {
    std::string name;
    LossType type;
};

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


Containers::Array<ComboElementFunctional> makeComboMapFunctionals() {
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




void drawDirichletEnergyOptions(DirichletEnergy& f){
    ImGui::Text("Dirichlet Energy");
}

void drawDoubleWellPotentialOptions(DoubleWellPotential& f){
    ImGui::Text("Double Well Potential");
}

void drawAreaRegularization1Options(AreaRegularizer1& f){
    constexpr double min = 0.f, max = 1.;
    ImGui::SliderScalar("Area Ratio", ImGuiDataType_Double, &f.areaRatio, &min, &max, "%.5f", 1.0f);
}

void drawAreaRegularization2Options(AreaRegularizer2& f){
    constexpr double min = 0.f, max = 1.;
    ImGui::SliderScalar("Area Ratio", ImGuiDataType_Double, &f.areaRatio, &min, &max, "%.5f", 1.0f);
}


