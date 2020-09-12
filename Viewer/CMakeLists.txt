
#find_package(MagnumPlugins REQUIRED PrimitiveImporter)

corrade_add_resource(Viewer_RESOURCES resources.conf)

add_library(Viewer STATIC
        Viewer.cpp
        Viewer.h
        ArcBall.cpp
        ArcBall.h
        KDTree.h
        ${viewer_RESOURCES})


target_include_directories(Viewer PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(Viewer
        PUBLIC
            Phasefield::Utilities
            ScopedTimer::ScopedTimer
            Phasefield::Cost
            Phasefield::Optimization
            Phasefield::Visualization
            Phasefield::ImGuiFS
            Phasefield::VisualizationProxy

            Magnum::GL
            Magnum::SceneGraph
            Magnum::Magnum
            Magnum::MeshTools
            Magnum::SceneGraph
            Magnum::Application
            Magnum::Shaders
            Magnum::DebugTools
            MagnumIntegration::ImGui
            ImGui::SourcesMiscCpp
        )

add_library(Phasefield::Viewer ALIAS Viewer)