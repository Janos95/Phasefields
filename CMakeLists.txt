cmake_minimum_required(VERSION 3.13)
project(phasefields)

set(CMAKE_CXX_STANDARD 17)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/modules)
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_CURRENT_LIST_DIR}/conan_modules)

find_package(Eigen3 CONFIG REQUIRED)
find_package(Ceres CONFIG REQUIRED)
find_package(Folly CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(OpenMesh CONFIG REQUIRED)
find_package(OpenMP CONFIG)
find_package(libigl CONFIG REQUIRED)
find_package(Corrade CONFIG REQUIRED PluginManager Containers TestSuite)
find_package(Magnum CONFIG REQUIRED GL SceneGraph MeshTools Sdl2Application Shaders Primitives)
find_package(MagnumPlugins CONFIG REQUIRED AssimpImporter)
find_package(MagnumIntegration CONFIG REQUIRED ImGui Eigen)
find_package(Boost CONFIG REQUIRED)
#find_package(ParallelSTL CONFIG REQUIRED)
#find_package(Enoki CONFIG REQUIRED)

#set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
#list(APPEND CMAKE_CXX_FLAGS "-fsanitize=address -fsanitize=undefined")
list(APPEND CMAKE_CXX_FLAGS "-ffast-math")

add_subdirectory(scoped_timer)
add_subdirectory(viewer)
add_subdirectory(cost)

find_package(PkgConfig REQUIRED)
pkg_check_modules(Ipopt REQUIRED IMPORTED_TARGET ipopt)

#set(CMAKE_CXX_STANDARD_LIBRARIES tbb)
#set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES /usr/local/include/c++/10.0.1)

add_executable(segmentation
        main.cpp
        color_callback.hpp
        my_configurator.hpp
        color_callback.cpp
        optimization_context.cpp
        optimization_context.hpp
        colormaps.hpp
        colormaps.cpp)

target_link_libraries(segmentation PRIVATE
        PkgConfig::Ipopt
        phasefields::viewer
        phasefields::scoped_timer
        phasefields::cost
        MagnumIntegration::Eigen
        Magnum::Primitives
        fmt::fmt
        ceres
        MagnumPlugins::AssimpImporter
        )

        
target_compile_definitions(segmentation PRIVATE CONF_PATH="${CMAKE_CURRENT_LIST_DIR}/phasefields.conf")

add_subdirectory(tests)
