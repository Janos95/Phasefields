

add_library(io
        load_mesh_data.cpp
        load_mesh_data.hpp
        write_mesh.cpp
        write_mesh.hpp)


target_include_directories(io PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(io
        PUBLIC
            Magnum::Trade
            Eigen3::Eigen
        PRIVATE
            phasefields::scoped_timer
            Corrade::Containers
            Corrade::PluginManager
            Magnum::AnySceneImporter
        )

add_library(phasefields::io ALIAS io)