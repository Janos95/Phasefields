

add_library(viewer
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
        draw_callbacks.hpp)

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
        PRIVATE
            Boost::Boost
        )


add_library(phasefields::viewer ALIAS viewer)