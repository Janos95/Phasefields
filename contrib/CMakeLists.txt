# Corrade
set(BUILD_PLUGINS_STATIC ON CACHE BOOL "" FORCE)
set(BUILD_STATIC ON CACHE BOOL "" FORCE)
set(WITH_CONTAINERS ON CACHE BOOL "" FORCE)
set(WITH_UTILITY ON CACHE BOOL "" FORCE)

add_subdirectory(corrade EXCLUDE_FROM_ALL)
find_package(Corrade) # bring in CORRADE_TARGET_* macros

# Magnum + Magnum-Plugins (Assimp)
if (PHASEFIELD_WITH_IO)
    set(WITH_ANYIMAGEIMPORTER ON CACHE BOOL "" FORCE)
    set(WITH_ASSIMPIMPORTER ON CACHE BOOL "" FORCE)
endif ()

set(WITH_STANFORDSCENECONVERTER ON CACHE BOOL "" FORCE)
set(WITH_STANFORDIMPORTER ON CACHE BOOL "" FORCE)

if (CORRADE_TARGET_EMSCRIPTEN)
    set(WITH_EMSCRIPTENAPPLICATION ON CACHE BOOL "" FORCE)
else()
    # GLFW
    set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
    set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)
    add_subdirectory(glfw EXCLUDE_FROM_ALL)

    set(WITH_GLFWAPPLICATION ON CACHE BOOL "" FORCE)
endif ()

add_subdirectory(magnum EXCLUDE_FROM_ALL)
add_subdirectory(magnum-plugins EXCLUDE_FROM_ALL)

# Magnum-Integration
set(IMGUI_DIR "${CMAKE_CURRENT_SOURCE_DIR}/imgui")
set(WITH_IMGUI ON CACHE BOOL "" FORCE)
add_subdirectory(magnum-integration EXCLUDE_FROM_ALL)
find_package(MagnumIntegration REQUIRED ImGui)

# ImPlot
add_library(implot implot/implot.h implot/implot_internal.h implot/implot.cpp implot/implot_items.cpp)
target_include_directories(implot PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/implot)
target_link_libraries(implot PUBLIC MagnumIntegration::ImGui)
add_library(ImPlot::ImPlot ALIAS implot)

# lbfgs
add_subdirectory(liblbfgs EXCLUDE_FROM_ALL)

# Bounding Volume Hierarchy
add_subdirectory(bvh EXCLUDE_FROM_ALL)

# Eigen
add_library(Eigen INTERFACE)
target_include_directories(Eigen INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/eigen/)
add_library(Eigen3::Eigen ALIAS Eigen)

set(WITH_DEMOS OFF CACHE BOOL "" FORCE)
add_subdirectory(SuiteSparse)

