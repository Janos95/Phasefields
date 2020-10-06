set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

add_library(Optimization STATIC
        Functional.h
        Functional.hpp
        Functional.cpp
        LossFunction.cpp
        LossFunctions.h
        Solver.cpp
        Solver.h
        Optimization.h
        Tag.h
        Tag.cpp
        RecursiveProblem.h
        RecursiveProblem.cpp
        Tree.h
        Tree.cpp
        )

target_include_directories(Optimization PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_link_libraries(Optimization PRIVATE
        Ceres::ceres
        ScopedTimer::ScopedTimer
        ImGui::ImGui
        Phasefield::Viewer
        )

target_link_libraries(Optimization PUBLIC
        Phasefield::Mesh
        Phasefield::VisualizationProxy
        Corrade::Containers
        )

add_library(Phasefield::Optimization ALIAS Optimization)

if (NOT CORRADE_TARGET_EMSCRIPTEN)
    target_link_libraries(Optimization PRIVATE
            TBB::tbb
            OpenMP::OpenMP_CXX
            )

    target_link_libraries(Optimization PUBLIC
            PkgConfig::Ipopt
            PkgConfig::Adolc
            )
endif ()


