//
// Created by janos on 26.05.20.
//

#include "meta_data.hpp"
#include "normalizeInto.hpp"

#include <Magnum/Trade/MeshData.h>
#include <imgui.h>
#include <mutex>

void MetaData::drawImGuiOptions(bool&, DrawableType&, bool&){
    ImGui::Text(name.c_str());
}

GradientMetaData::GradientMetaData(Mg::Trade::MeshData& md, VisualizationFlags& u, std::mutex& m) :
    meshData(md), update(u), mutex(m)
{
}

void GradientMetaData::visualizeGradient(Cr::Containers::ArrayView<const Mg::Double> const& gradient) {

    std::lock_guard l(mutex);
    if(flags & VisualizationFlag::Gradient){
        auto coords = meshData.mutableAttribute(Mg::Trade::MeshAttribute::TextureCoordinates);
        auto xcoords = Cr::Containers::arrayCast<2, Mg::Float>(coords).template slice<1>();
        normalizeInto(gradient, xcoords);
        update = VisualizationFlag::Gradient;
    }
}

void GradientMetaData::drawImGuiOptions(bool& makeExclusive, DrawableType& type, bool& evaluateProblem){
    ImGui::Text(this->name.c_str());
    bool drawGrad = bool(flags & VisualizationFlag::Gradient);
    if(ImGui::Checkbox("Gradient", &drawGrad) ){
        std::lock_guard l(mutex);
        if(drawGrad){
            type = DrawableType::PhongDiffuse;
            flags = VisualizationFlag::Gradient;
            makeExclusive = true;
            evaluateProblem = true;
        } else {
            flags = {};
        }
    }
}

void AreaMetaData::drawImGuiOptions(bool& makeExclusive, DrawableType& type, bool& evaluateProblem) {
    GradientMetaData::drawImGuiOptions(makeExclusive,type,evaluateProblem);
    constexpr Mg::Double min = 0., max = 1.;
    ImGui::SliderScalar("Area Ratio", ImGuiDataType_Double, &areaRatio, &min, &max, "%.5f", 1.0f);
}