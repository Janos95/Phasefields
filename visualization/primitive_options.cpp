//
// Created by janos on 18.05.20.
//

#include "primitive_options.hpp"
#include <imgui.h>

bool PolygonizationOptions::drawImGui() {
    constexpr static float minAngle = 10.f, maxAngle = 60.f;
    constexpr static float minBoundingRadius = 1., maxBoundingRadius = 10.f;
    constexpr static float minDistanceBound = 0.001f, maxDistanceBound = .5f;
    constexpr static float minRadiusBound = 0.001f, maxRadiusBound = .5f;
    bool refresh = false;
    refresh |= ImGui::SliderScalar("Angle Bound", ImGuiDataType_Float, &angleBound, &minAngle, &maxAngle, "%f");
    refresh |= ImGui::SliderScalar("Bounding Sphere Radius", ImGuiDataType_Float, &boundingSphereRadius, &minBoundingRadius, &maxBoundingRadius, "%f");
    refresh |= ImGui::SliderScalar("Distance Bound", ImGuiDataType_Float, &distanceBound, &minDistanceBound, &maxDistanceBound, "%.2e");
    refresh |= ImGui::SliderScalar("Radius Bound", ImGuiDataType_Float, &radiusBound, &minRadiusBound, &maxRadiusBound, "%.2e");
    return refresh;
}

