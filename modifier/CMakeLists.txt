find_package(Corrade CONFIG REQUIRED PluginManager Containers TestSuite)
find_package(Magnum CONFIG REQUIRED GL SceneGraph MeshTools Sdl2Application Shaders Primitives DebugTools)
find_package(MagnumPlugins CONFIG REQUIRED AssimpImporter)
find_package(MagnumIntegration CONFIG REQUIRED ImGui Eigen)
find_package(ImGui REQUIRED COMPONENTS SourcesMiscCpp)
find_package(Eigen3 CONFIG REQUIRED)

add_library(modifier STATIC
        optimization_context.cpp
        optimization_context.hpp
        brush.cpp
        brush.hpp
        subdivision.cpp
        subdivision.hpp
        mesh_io.cpp
        mesh_io.hpp
        primitives.cpp
        primitives.hpp
        shader_options.hpp
        shader_options.cpp
        custom_widgets.cpp
        custom_widgets.hpp
        phasefield_data.hpp
        )

target_include_directories(modifier PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_link_libraries(modifier PUBLIC
        phasefields::cost
        )
target_link_libraries(modifier PRIVATE
        phasefields::scoped_timer
        phasefields::exact_geodesics
        phasefields::polygonizer
        Eigen3::Eigen
        Magnum::Primitives
        MagnumPlugins::AssimpImporter
        ImGui::SourcesMiscCpp
        )

add_library(phasefields::modifier ALIAS modifier)
