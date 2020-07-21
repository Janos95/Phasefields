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


add_library(visualization_proxy STATIC
        visualization_proxy.cpp
        visualization_proxy.hpp
        )

target_include_directories(visualization_proxy PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(visualization_proxy
        PRIVATE
            Magnum::GL
            Magnum::Magnum
            phasefields::viewer
            PUBLIC
            phasefields::utilities
        )

set_property(TARGET visualization_proxy PROPERTY POSITION_INDEPENDENT_CODE ON)
add_library(phasefields::visualization_proxy ALIAS visualization)
