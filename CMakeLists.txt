cmake_minimum_required(VERSION 3.13)
project(Phasefield)

set(CMAKE_CXX_STANDARD 20)

#add_subdirectory(contrib)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/modules)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#find_package(TBB CONFIG REQUIRED)
#find_package(ipopt REQUIRED)

#set_compile_options("-Wno-old-style-cast")
#list(INSERT CMAKE_CXX_FLAGS 0 "-Wno-old-style-cast")
#set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)


find_package(Corrade REQUIRED Containers Utility)
find_package(Magnum REQUIRED GL MeshTools Shaders Primitives DebugTools)
find_package(MagnumPlugins REQUIRED PrimitiveImporter)
find_package(MagnumIntegration REQUIRED ImGui)
find_package(Eigen3)
find_package(Ceres)

if (NOT CORRADE_TARGET_EMSCRIPTEN)
    find_package(TBB CONFIG REQUIRED)
    find_package(PkgConfig REQUIRED)
    find_package(OpenMP)
    find_package(SuiteSparse REQUIRED UMFPACK)
    find_package(MagnumPlugins REQUIRED AssimpImporter)
    find_package(Magnum REQUIRED Sdl2Application)

    pkg_check_modules(Ipopt REQUIRED IMPORTED_TARGET ipopt)
    pkg_check_modules(Adolc REQUIRED IMPORTED_TARGET adolc)

    pkg_check_modules(libavcodec REQUIRED IMPORTED_TARGET libavcodec)
    pkg_check_modules(libavutil REQUIRED IMPORTED_TARGET libavutil)
    pkg_check_modules(libswscale REQUIRED IMPORTED_TARGET libswscale)
else()
    find_package(Magnum REQUIRED EmscriptenApplication)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_subdirectory(Utilities)
add_subdirectory(ScopedTimer)
add_subdirectory(Visualization)
#add_subdirectory(polygonizer)
#add_subdirectory(math)
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
