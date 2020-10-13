corrade_add_resource(Viewer_RESOURCES resources.conf)
add_library(Viewer STATIC
        Viewer.cpp
        Viewer.h
        ArcBall.cpp
        ArcBall.h
        KDTree.h
        KDTree.cpp
        DepthReinterpretShader.cpp
        DepthReinterpretShader.h
        ${viewer_RESOURCES})

target_include_directories(Viewer PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_link_libraries(Viewer
        PUBLIC
            Phasefield::Utilities
            ScopedTimer::ScopedTimer
            Phasefield::Cost
            Phasefield::Optimization
            Phasefield::Visualization
            Phasefield::VisualizationProxy

            Magnum::GL
            Magnum::Magnum
            Magnum::MeshTools
            Magnum::Application
            Magnum::Shaders
            Magnum::DebugTools
            MagnumPlugins::PrimitiveImporter
            MagnumIntegration::ImGui
            ImPlot::ImPlot
            FastBVH::FastBVH
        )

add_library(Phasefield::Viewer ALIAS Viewer)

if (PHASEFIELD_WITH_ASSIMP)
    target_link_libraries(Viewer PUBLIC  MagnumPlugins::AssimpImporter)
    target_compile_definitions(Viewer PUBLIC PHASEFIELD_WITH_ASSIMP)
endif ()

if (PHASEFIELD_WITH_FFMPEG)
    target_link_libraries(Viewer PUBLIC Phasefield::VideoSaver)
endif()
