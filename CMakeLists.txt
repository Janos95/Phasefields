cmake_minimum_required(VERSION 3.13)
project(phasefields)

set(CMAKE_CXX_STANDARD 20)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/modules)

find_package(Eigen3 CONFIG REQUIRED)
find_package(Ceres CONFIG REQUIRED)
find_package(Folly CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(goast CONFIG REQUIRED)
find_package(OpenMesh CONFIG REQUIRED)
find_package(OpenMP CONFIG)
find_package(libigl CONFIG REQUIRED)
find_package(Corrade CONFIG REQUIRED PluginManager Containers TestSuite)
find_package(Magnum CONFIG REQUIRED GL SceneGraph MeshTools Sdl2Application Shaders AnySceneImporter Primitives)
find_package(MagnumIntegration CONFIG REQUIRED ImGui Eigen)
#find_package(Enoki CONFIG REQUIRED)

#set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fsanitize=undefined")

add_subdirectory(scoped_timer)
add_subdirectory(viewer)
add_subdirectory(io)
add_subdirectory(cost)

find_package(PkgConfig REQUIRED)
pkg_check_modules(Ipopt REQUIRED IMPORTED_TARGET ipopt)

add_executable(segmentation
        main.cpp
        color_callback.hpp
        phasefield_initialization.hpp
        my_configurator.hpp
        color_callback.cpp
        colormaps.hpp
        colormaps.cpp)

target_link_libraries(segmentation PRIVATE
        goast::goast
        PkgConfig::Ipopt
        phasefields::viewer
        phasefields::scoped_timer
        phasefields::cost
        phasefields::io
        MagnumIntegration::Eigen
        Magnum::Primitives
        OpenMeshCore
        ceres)


add_subdirectory(tests)
