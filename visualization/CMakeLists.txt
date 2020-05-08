add_library(visualization STATIC
        paths.hpp
        paths.cpp
        types.hpp
        upload.cpp
        upload.hpp
        drawables.cpp
        drawables.hpp
        primitive_options.hpp
        )
target_include_directories(visualization PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(visualization
        PUBLIC
        Magnum::GL
        Magnum::SceneGraph
        Magnum::Magnum
        Magnum::MeshTools
        Magnum::SceneGraph
        )
set_property(TARGET visualization PROPERTY POSITION_INDEPENDENT_CODE ON)
add_library(phasefields::visualization ALIAS visualization)
