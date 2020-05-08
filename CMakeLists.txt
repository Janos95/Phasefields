cmake_minimum_required(VERSION 3.13)
project(phasefields)

set(CMAKE_CXX_STANDARD 20)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/modules)
#list(APPEND CMAKE_PREFIX_PATH ${CMAKE_CURRENT_LIST_DIR}/conan_modules)

find_package(Eigen3 CONFIG REQUIRED)
find_package(Corrade CONFIG REQUIRED PluginManager Containers TestSuite)
find_package(Magnum CONFIG REQUIRED GL SceneGraph MeshTools Sdl2Application Shaders Primitives DebugTools)
find_package(MagnumPlugins CONFIG REQUIRED AssimpImporter)
find_package(MagnumIntegration CONFIG REQUIRED ImGui Eigen)
find_package(ImGui REQUIRED COMPONENTS SourcesMiscCpp)
find_package(TBB CONFIG REQUIRED)
#find_package(ParallelSTL CONFIG REQUIRED)
#find_package(Enoki CONFIG REQUIRED)

#set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
#list(APPEND CMAKE_CXX_FLAGS "-fsanitize=address -fsanitize=undefined")
list(APPEND CMAKE_CXX_FLAGS "-frounding-math")

add_subdirectory(utilities)
add_subdirectory(scoped_timer)
add_subdirectory(visualization)
add_subdirectory(polygonizer)
add_subdirectory(cost)
add_subdirectory(exact_geodesics)
add_subdirectory(shaders)
add_subdirectory(viewer)


#set(CMAKE_CXX_STANDARD_LIBRARIES tbb)
#set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES /usr/local/include/c++/10.0.1)

#add_executable(segmentation
#        main.cpp
#        )

#target_link_libraries(segmentation PRIVATE
#        phasefields::viewer
#        )
#add_subdirectory(tests)
