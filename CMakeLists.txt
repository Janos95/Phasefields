cmake_minimum_required(VERSION 3.13)
project(Phasefield)

set(CMAKE_CXX_STANDARD 20)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/modules)
#list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/conan_modules)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_FLAGS_RELEASE "-ffast-math")

find_package(Eigen3 CONFIG REQUIRED)
find_package(Corrade CONFIG REQUIRED PluginManager Containers TestSuite)
find_package(Magnum CONFIG REQUIRED GL SceneGraph MeshTools Sdl2Application Shaders Primitives DebugTools)
find_package(MagnumPlugins CONFIG REQUIRED AssimpImporter)
find_package(MagnumIntegration CONFIG REQUIRED ImGui Eigen)
find_package(ImGui REQUIRED COMPONENTS SourcesMiscCpp)
find_package(TBB CONFIG REQUIRED)
#find_package(ipopt REQUIRED)

#set_compile_options("-Wno-old-style-cast")
#list(INSERT CMAKE_CXX_FLAGS 0 "-Wno-old-style-cast")
#set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_subdirectory(Utilities)
add_subdirectory(ImGuiFS)
add_subdirectory(ScopedTimer)
add_subdirectory(Visualization)
#add_subdirectory(polygonizer)
#add_subdirectory(math)
#add_subdirectory(Shaders)
add_subdirectory(Optimization)
add_subdirectory(Mesh)
add_subdirectory(Cost)
add_subdirectory(VideoSaver)
add_subdirectory(Viewer)

#set(CMAKE_CXX_STANDARD_LIBRARIES tbb)
#set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES /usr/local/include/c++/10.0.1)

add_executable(Phasefield main.cpp)

target_link_libraries(Phasefield PRIVATE Phasefield::Viewer)
#set_property(TARGET Phasefield PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
#add_subdirectory(tests)
