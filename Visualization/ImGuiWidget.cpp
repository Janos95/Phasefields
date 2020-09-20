//
// Created by janos on 08.04.20.
// taken from some imgui github issue
//

#include "ImGuiWidgets.h"

#include <imgui.h>
#include <imgui_internal.h>

bool toggleButton(const char* str_id, bool* v)
{

    bool clicked = false;
    ImVec2 p = ImGui::GetCursorScreenPos();
    ImDrawList* draw_list = ImGui::GetWindowDrawList();

    float height = ImGui::GetFrameHeight();
    float width = height * 1.55f;
    float radius = height * 0.50f;

    ImGui::InvisibleButton(str_id, ImVec2(width, height));
    if (ImGui::IsItemClicked()){
        clicked = true;
        *v = !*v;
    }

    float t = *v ? 1.0f : 0.0f;

    ImGuiContext& g = *GImGui;
    float ANIM_SPEED = 0.08f;
    if (g.LastActiveId == g.CurrentWindow->GetID(str_id))// && g.LastActiveIdTimer < ANIM_SPEED)
    {
        float t_anim = ImSaturate(g.LastActiveIdTimer / ANIM_SPEED);
        t = *v ? (t_anim) : (1.0f - t_anim);
    }

    ImU32 col_bg;
    if (ImGui::IsItemHovered())
        col_bg = ImGui::GetColorU32(ImLerp(ImVec4(0.78f, 0.78f, 0.78f, 1.0f), ImVec4(0.64f, 0.83f, 0.34f, 1.0f), t));
    else
        col_bg = ImGui::GetColorU32(ImLerp(ImVec4(0.85f, 0.85f, 0.85f, 1.0f), ImVec4(0.56f, 0.83f, 0.26f, 1.0f), t));

    draw_list->AddRectFilled(p, ImVec2(p.x + width, p.y + height), col_bg, height * 0.5f);
    draw_list->AddCircleFilled(ImVec2(p.x + radius + t * (width - radius * 2.0f), p.y + radius), radius - 1.5f, IM_COL32(255, 255, 255, 255));

    return clicked;
}

bool dragDoubleRange2(const char* label, double* v_current_min, double* v_current_max, float v_speed, double v_min, double v_max, const char* format, const char* format_max, float power)
{
    ImGuiWindow* window = ImGui::GetCurrentWindow();
    if (window->SkipItems)
        return false;

    ImGuiContext& g = *GImGui;
    ImGui::PushID(label);
    ImGui::BeginGroup();
    ImGui::PushMultiItemsWidths(2, ImGui::CalcItemWidth());

    auto min = (v_min >= v_max) ? -DBL_MAX : v_min;
    auto max = (v_min >= v_max) ? *v_current_max : ImMin(v_max, *v_current_max);
    bool value_changed = ImGui::DragScalar("##min", ImGuiDataType_Double, v_current_min, v_speed, &min, &max, format, power);
    ImGui::PopItemWidth();
    ImGui::SameLine(0, g.Style.ItemInnerSpacing.x);
    min = (v_min >= v_max) ? *v_current_min : ImMax(v_min, *v_current_min);
    max = (v_min >= v_max) ? FLT_MAX : v_max;
    value_changed |= ImGui::DragScalar("##max", ImGuiDataType_Double, v_current_max, v_speed, &min, &max, format_max ? format_max : format, power);
    ImGui::PopItemWidth();
    ImGui::SameLine(0, g.Style.ItemInnerSpacing.x);

    ImGui::TextEx(label, ImGui::FindRenderedTextEnd(label));
    ImGui::EndGroup();
    ImGui::PopID();
    return value_changed;
}

//ImRect RenderTree(Node* n)
//{
//    const bool recurse = ImGui::TreeNode(...);
//    const ImRect nodeRect = ImRect(ImGui::GetItemRectMin(), ImGui::GetItemRectMax());
//
//    if (recurse)
//    {
//        const ImColor TreeLineColor = ImGui::GetColorU32(ImGuiCol_Text);
//        const float SmallOffsetX = 11.0f; //for now, a hardcoded value; should take into account tree indent size
//        ImDrawList* drawList = ImGui::GetWindowDrawList();
//
//        ImVec2 verticalLineStart = ImGui::GetCursorScreenPos();
//        verticalLineStart.x += SmallOffsetX; //to nicely line up with the arrow symbol
//        ImVec2 verticalLineEnd = verticalLineStart;
//
//        for (Node* child : *n)
//        {
//            const float HorizontalTreeLineSize = 8.0f; //chosen arbitrarily
//            const ImRect childRect = RenderTree(child);
//            const float midpoint = (childRect.Min.y + childRect.Max.y) / 2.0f;
//            drawList->AddLine(ImVec2(verticalLineStart.x, midpoint), ImVec(verticalLineStart.x + HorizontalTreeLineSize, midpoint), TreeLineColor);
//            verticalLineEnd.y = midpoint;
//        }
//
//        drawList->AddLine(verticalLineStart, verticalLineEnd, TreeLineColor);
//    }
//
//    return nodeRect;
//}
