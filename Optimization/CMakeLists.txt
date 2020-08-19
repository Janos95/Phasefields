find_package(Ceres CONFIG REQUIRED)
find_package(TBB CONFIG REQUIRED)

find_package(PkgConfig REQUIRED)

list(APPEND CMAKE_PREFIX_PATH /home/janos/libraries/prefix/)
set(PKG_CONFIG_USE_CMAKE_PREFIX_PATH ON)

pkg_check_modules(Adolc REQUIRED IMPORTED_TARGET adolc)

set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

add_library(Optimization STATIC
        Functional.h
        Functional.cpp
        LossFunction.cpp
        LossFunctions.h
        Solver.cpp
        Solver.h
        Optimization.h
        Functional.hpp
        Tag.h
        Tag.cpp
        RecursiveProblem.h
        RecursiveProblem.cpp
        Tree.h
        Tree.cpp)

target_include_directories(Optimization PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_link_libraries(Optimization PRIVATE
        Ceres::ceres
        ipopt::ipopt
        TBB::tbb
        ScopedTimer::ScopedTimer
        ImGui::ImGui
        )

target_link_libraries(Optimization PUBLIC
        PkgConfig::Adolc
        Phasefield::VisualizationProxy
        Corrade::Containers
        )

set_property(TARGET Optimization PROPERTY POSITION_INDEPENDENT_CODE ON)

add_library(Phasefield::Optimization ALIAS Optimization)


