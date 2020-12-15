cmake_minimum_required(VERSION 3.13)
project(Phasefield)

cmake_policy(SET CMP0077 NEW)

set(MKL_ROOT /opt/intel/mkl)
set(MKL_PARALLEL ON)
set(MKL_64BIT ON)

set(CMAKE_CXX_STANDARD 20)
#set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
#set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/modules)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
#set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

option(PHASEFIELD_WITH_ADOLC "Enable automatica gradient testing" ON)
option(PHASEFIELD_WITH_CERES "Enables ceres for optimization" ON)
option(PHASEFIELD_WITH_IPOPT "Enable ipopt for consrained optimization" ON)
option(PHASEFIELD_WITH_TBB "Enable threading using intel tbb" ON)
option(PHASEFIELD_WITH_SUITESPARSE "Use umfpack from suite sparse for faster decompositions" ON)
option(PHASEFIELD_WITH_IO "Import with assimp, export with vtk to vts (paraview), ply and videos using ffmpeg" ON)

add_subdirectory(contrib)

find_package(MKL)
find_package(PkgConfig REQUIRED)
find_package(Corrade)

if(PHASEFIELD_WITH_IO)
    find_package(MagnumPlugins REQUIRED AssimpImporter)
    find_package(MagnumPlugins REQUIRED StanfordSceneConverter)
    find_package(VTK REQUIRED CommonCore CommonDataModel IOXML)

    pkg_check_modules(libavcodec REQUIRED IMPORTED_TARGET libavcodec)
    pkg_check_modules(libavutil REQUIRED IMPORTED_TARGET libavutil)
    pkg_check_modules(libswscale REQUIRED IMPORTED_TARGET libswscale)
    add_subdirectory(VideoSaver)
endif()

if(PHASEFIELD_WITH_ADOLC)
    #message(FATAL_ERROR "finding adolc")
    pkg_check_modules(Adolc REQUIRED IMPORTED_TARGET adolc)
endif()

if(PHASEFIELD_WITH_IPOPT)
    pkg_check_modules(Ipopt REQUIRED IMPORTED_TARGET ipopt)
endif()

if(PHASEFIELD_WITH_CERES)
    find_package(Ceres)
endif()

if(PHASEFIELD_WITH_TBB)
    find_package(TBB)
endif()

if(PHASEFIELD_WITH_SUITESPARSE)
    find_package(SuiteSparse REQUIRED UMFPACK)
endif()

if (CORRADE_TARGET_EMSCRIPTEN)
    find_package(Magnum REQUIRED EmscriptenApplication)
    #if(NOT TARGET OpenGLES2::OpenGLES2)
    #    message(FATAL_ERROR "gles2 not found")
    #elseif()
    #    message(FATAL_ERROR "gles2 found")
    #endif()
else()
    find_package(Magnum REQUIRED Sdl2Application)
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
add_subdirectory(Viewer)

add_executable(Phasefield main.cpp)

target_link_libraries(Phasefield PRIVATE Phasefield::Viewer)
#set_property(TARGET Phasefield PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
#add_subdirectory(tests)

# For Emscripten copy the boilerplate next to the executable so it can be run
# directly from the build dir; provide an install target as well
if(CORRADE_TARGET_EMSCRIPTEN)
    add_custom_command(TARGET Phasefield POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${MAGNUM_EMSCRIPTENAPPLICATION_JS}
            ${MAGNUM_WEBAPPLICATION_CSS}
            $<TARGET_FILE_DIR:Phasefield>
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_CURRENT_SOURCE_DIR}/Phasefield.html
            $<TARGET_FILE_DIR:Phasefield>/Phasefield.html)

    install(FILES Phasefield.html DESTINATION ${MAGNUM_DEPLOY_PREFIX}/Phasefield RENAME index.html)
    install(TARGETS Phasefield DESTINATION ${MAGNUM_DEPLOY_PREFIX}/Phasefield)
    install(FILES
            ${MAGNUM_EMSCRIPTENAPPLICATION_JS}
            ${MAGNUM_WEBAPPLICATION_CSS}
            DESTINATION ${MAGNUM_DEPLOY_PREFIX}/Phasefield)
    install(FILES
            $<TARGET_FILE_DIR:Phasefield>/Phasefield.js.mem
            $<TARGET_FILE_DIR:Phasefield>/Phasefield.wasm
            DESTINATION ${MAGNUM_DEPLOY_PREFIX}/Phasefield OPTIONAL)
endif()