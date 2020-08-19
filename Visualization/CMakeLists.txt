add_library(Visualization STATIC
        Paths.h
        Paths.cpp
        Upload.cpp
        Upload.h
        primitive_options.hpp
        primitive_options.cpp
        custom_widgets.hpp
        custom_widgets.cpp
        )
target_include_directories(Visualization PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(Visualization
        PUBLIC
            Magnum::GL
            Magnum::SceneGraph
            Magnum::Magnum
            Magnum::MeshTools
            Magnum::SceneGraph
            MagnumIntegration::ImGui
            ImGui::SourcesMiscCpp
            Phasefield::Utilities
        )
set_property(TARGET Visualization PROPERTY POSITION_INDEPENDENT_CODE ON)
add_library(Phasefield::Visualization ALIAS Visualization)


add_library(VisualizationProxy STATIC
        VisualizationProxy.cpp
        VisualizationProxy.h
        )

target_include_directories(VisualizationProxy PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(VisualizationProxy
        PRIVATE
            Magnum::GL
            Magnum::Magnum
            Phasefield::Viewer
        PUBLIC
            Phasefield::Utilities
        )

set_property(TARGET VisualizationProxy PROPERTY POSITION_INDEPENDENT_CODE ON)
add_library(Phasefield::VisualizationProxy ALIAS VisualizationProxy)
