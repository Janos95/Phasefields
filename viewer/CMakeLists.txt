

corrade_add_resource(viewer_RESOURCES resources.conf)

add_library(viewer SHARED
        viewer.cpp
        viewer.hpp
        arc_ball.cpp
        arc_ball.hpp
        arc_ball_camera.hpp
        camera.hpp
        upload.cpp
        upload.hpp
        drawables.cpp
        drawables.hpp
        primitives.hpp
        primitives.cpp
        shader_options.hpp
        shader_options.cpp
        subdivision.hpp
        subdivision.cpp
        ${viewer_RESOURCES})

target_include_directories(viewer PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(viewer
        PUBLIC
            phasefields::utilities
            phasefields::scoped_timer
            phasefields::cost
            phasefields::exact_geodesics
            phasefields::polygonizer
            Magnum::GL
            Magnum::SceneGraph
            Magnum::Magnum
            Magnum::MeshTools
            Magnum::SceneGraph
            Magnum::Application
            Magnum::Shaders
            Magnum::DebugTools
            MagnumIntegration::ImGui
        )

add_library(phasefields::viewer ALIAS viewer)