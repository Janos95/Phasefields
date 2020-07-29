find_package(Ceres CONFIG REQUIRED)
find_package(TBB CONFIG REQUIRED)

find_package(PkgConfig REQUIRED)

list(APPEND CMAKE_PREFIX_PATH /home/janos/libraries/prefix/)
set(PKG_CONFIG_USE_CMAKE_PREFIX_PATH ON)

pkg_check_modules(Adolc REQUIRED IMPORTED_TARGET adolc)

set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

add_library(optimization STATIC
        Functional.h
        Functional.cpp
        LossFunction.cpp
        LossFunctions.h
        Solver.cpp
        Solver.h
        SparseMatrix.h
        Optimization.h
        Functional.hpp
        SparseMatrix.cpp
        tag.h
        tag.cpp
        RecursiveProblem.h
        RecursiveProblem.cpp PhasefieldTree.cpp)

target_include_directories(optimization PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_link_libraries(optimization PRIVATE
        ceres
        ipopt::ipopt
        TBB::tbb
        phasefields::scoped_timer
        ImGui::ImGui
        )

target_link_libraries(optimization PUBLIC
        PkgConfig::Adolc
        phasefields::visualization_proxy
        Corrade::Containers
        )

set_property(TARGET optimization PROPERTY POSITION_INDEPENDENT_CODE ON)

add_library(phasefields::optimization ALIAS optimization)


