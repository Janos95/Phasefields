corrade_add_resource(Viewer_Rsc resources.conf)
corrade_add_resource(Experiments_Rsc data/resources.conf)
add_library(Viewer STATIC
        Viewer.cpp
        Viewer.h
        ArcBall.cpp
        ArcBall.h
        KDTree.h
        KDTree.cpp
        Bvh.cpp
        Bvh.h
        ${Viewer_Rsc}
        ${Experiments_Rsc}
        )

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
            MagnumIntegration::ImGui
            MagnumPlugins::StanfordImporter
            ImPlot::ImPlot
            bvh
        )

add_library(Phasefield::Viewer ALIAS Viewer)

if (PHASEFIELD_WITH_IO)
    target_link_libraries(Viewer PUBLIC
            MagnumPlugins::AssimpImporter
            MagnumPlugins::StanfordSceneConverter
            VTK::CommonCore
            VTK::CommonDataModel
            VTK::IOXML
            Phasefield::VideoSaver)
    target_compile_definitions(Viewer PUBLIC PHASEFIELD_WITH_IO)
endif ()

