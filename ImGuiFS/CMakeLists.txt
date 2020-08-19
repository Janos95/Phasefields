
add_library(ImGuiFS imguifilesystem.cpp imguifilesystem.h dirent_portable.h)
target_include_directories(ImGuiFS SYSTEM PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_compile_options(ImGuiFS PRIVATE -w)
target_link_libraries(ImGuiFS PUBLIC ImGui::ImGui)
add_library(Phasefield::ImGuiFS ALIAS ImGuiFS)