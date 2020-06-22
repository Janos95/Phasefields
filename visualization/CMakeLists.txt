add_library(visualization STATIC
        paths.hpp
        paths.cpp
        upload.cpp
        upload.hpp
        drawables.cpp
        drawables.hpp
        primitive_options.hpp
        primitive_options.cpp
        custom_widgets.hpp
        custom_widgets.cpp
        visualization_proxy.hpp
        visualization_proxy.cpp
        )
target_include_directories(visualization PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(visualization
        PUBLIC
        Magnum::GL
        Magnum::SceneGraph
        Magnum::Magnum
        Magnum::MeshTools
        Magnum::SceneGraph
        MagnumIntegration::ImGui
        ImGui::SourcesMiscCpp
        phasefields::utilities
        )
set_property(TARGET visualization PROPERTY POSITION_INDEPENDENT_CODE ON)
add_library(phasefields::visualization ALIAS visualization)
