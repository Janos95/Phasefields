
set(BUILD_PLUGINS_STATIC ON CACHE BOOL "" FORCE)
set(BUILD_STATIC ON CACHE BOOL "" FORCE)

set(WITH_CONTAINERS ON CACHE BOOL "" FORCE)
set(WITH_UTILITY ON CACHE BOOL "" FORCE)

add_subdirectory(corrade EXCLUDE_FROM_ALL)
find_package(Corrade REQUIRED Containers Utility)

if (PHASEFIELD_WITH_IO)
    set(WITH_ANYIMAGEIMPORTER ON CACHE BOOL "" FORCE)
    set(WITH_ASSIMPIMPORTER ON CACHE BOOL "" FORCE)
endif ()

set(WITH_STANFORDSCENECONVERTER ON CACHE BOOL "" FORCE)
set(WITH_STANFORDIMPORTER ON CACHE BOOL "" FORCE)

if (CORRADE_TARGET_EMSCRIPTEN)
    set(WITH_EMSCRIPTENAPPLICATION ON CACHE BOOL "" FORCE)
else()
    set(WITH_SDL2APPLICATION ON CACHE BOOL "" FORCE)
endif ()

add_subdirectory(magnum EXCLUDE_FROM_ALL)
find_package(Magnum REQUIRED GL MeshTools Primitives Shaders Trade DebugTools)

add_subdirectory(magnum-plugins EXCLUDE_FROM_ALL)

set(IMGUI_DIR ${CMAKE_CURRENT_SOURCE_DIR}/imgui)
set(WITH_IMGUI ON CACHE BOOL "" FORCE)

add_subdirectory(magnum-integration EXCLUDE_FROM_ALL)
find_package(MagnumIntegration REQUIRED ImGui)

add_library(implot implot/implot.h implot/implot_internal.h implot/implot.cpp implot/implot_items.cpp)
target_include_directories(implot PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/implot)
target_link_libraries(implot PUBLIC MagnumIntegration::ImGui)
add_library(ImPlot::ImPlot ALIAS implot)

add_subdirectory(liblbfgs)

add_subdirectory(bvh)



