
add_library(imguifs imguifilesystem.cpp imguifilesystem.h dirent_portable.h)
target_include_directories(imguifs SYSTEM PUBLIC ${CMAKE_SOURCE_DIR})
target_compile_options(imguifs PRIVATE -w)
target_link_libraries(imguifs PUBLIC ImGui::ImGui)
add_library(phasefields::imguifs ALIAS imguifs)