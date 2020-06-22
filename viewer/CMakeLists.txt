

corrade_add_resource(viewer_RESOURCES resources.conf)

add_library(viewer SHARED
        viewer.cpp
        viewer.hpp
        arc_ball.cpp
        arc_ball.hpp
        arc_ball_camera.hpp
        primitives.hpp
        primitives.cpp
        subdivision.hpp
        subdivision.cpp
        optimization_context.hpp
        optimization_context.cpp
        imguifilesystem.cpp
        imguifilesystem.h
        ${viewer_RESOURCES})



target_include_directories(viewer PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(viewer
        PUBLIC
            phasefields::utilities
            phasefields::scoped_timer
            phasefields::cost
            phasefields::exact_geodesics
            phasefields::visualization
            #phasefields::polygonizer
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

add_library(phasefields::viewer ALIAS viewer)