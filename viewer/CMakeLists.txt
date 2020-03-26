

corrade_add_resource(viewer_RESOURCES resources.conf)

add_library(viewer SHARED
        viewer.cpp
        viewer.hpp
        scene_graph_node.cpp
        scene_graph_node.hpp
        arc_ball.cpp
        arc_ball.hpp
        arc_ball_camera.hpp
        scene.hpp
        scene.cpp
        object.hpp
        camera.hpp
        upload.cpp
        upload.hpp
        draw_callbacks.cpp
        draw_callbacks.hpp
        ${viewer_RESOURCES})

target_include_directories(viewer PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(viewer
        PUBLIC
            Magnum::GL
            Magnum::SceneGraph
            Magnum::Magnum
            Magnum::MeshTools
            Magnum::SceneGraph
            Magnum::Application
            Magnum::Shaders
            MagnumIntegration::ImGui
        )

add_library(phasefields::viewer ALIAS viewer)