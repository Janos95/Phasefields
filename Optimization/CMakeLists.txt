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
        LbfgsSolver.cpp
        LbfgsSolver.h
        )

target_include_directories(Optimization PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_link_libraries(Optimization PRIVATE
        LBFGS::lib
        ScopedTimer::ScopedTimer
        Phasefield::Viewer
        )

target_link_libraries(Optimization PUBLIC
        Phasefield::Mesh
        Phasefield::VisualizationProxy
        Corrade::Containers
        )

add_library(Phasefield::Optimization ALIAS Optimization)

if (PHASEFIELD_WITH_CERES)
    target_link_libraries(Optimization PRIVATE Ceres::ceres)
    target_compile_definitions(Optimization PUBLIC PHASEFIELD_WITH_CERES)
endif()

if (PHASEFIELD_WITH_IPOPT)
    target_link_libraries(Optimization PRIVATE PkgConfig::Ipopt)
    target_compile_definitions(Optimization PUBLIC PHASEFIELD_WITH_IPOPT)
endif()

if (PHASEFIELD_WITH_ADOLC)
    target_link_libraries(Optimization PRIVATE PkgConfig::Adolc)
    target_compile_definitions(Optimization PUBLIC PHASEFIELD_WITH_ADOLC)
endif()



