find_package(Ceres CONFIG REQUIRED)
find_package(TBB CONFIG REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(Ipopt REQUIRED IMPORTED_TARGET ipopt)
pkg_check_modules(Adolc REQUIRED IMPORTED_TARGET adolc)

set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

add_library(cost STATIC
        bfs.hpp
        c1_functions.hpp
        connectedness_constraint.hpp
        connectedness_constraint.h
        connectedness_constraint.cpp
        cost.h
        dijkstra.hpp
        fem.cpp
        fem.hpp
        functional.hpp
        graph_common.hpp
        loss_functions.cpp
        loss_functions.hpp
        meta_data.hpp
        meta_data.cpp
        modica_mortola.cpp
        modica_mortola.hpp
        problem.cpp
        problem.hpp
        small_array.hpp
        solver.cpp
        solver.hpp
        stopping_criteria.cpp
        stopping_criteria.hpp
        union_find.cpp
        union_find.hpp
        strang_rules.hpp
        phasefield_tree.hpp
        adolc_adapter.h
        adolc_adapter.hpp)

target_include_directories(cost PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_link_libraries(cost PRIVATE
        ceres
        TBB::tbb
        phasefields::scoped_timer
        )

target_link_libraries(cost PUBLIC
        PkgConfig::Ipopt
        PkgConfig::Adolc
        Eigen3::Eigen
        phasefields::utilities
        phasefields::visualization
        )

set_property(TARGET cost PROPERTY POSITION_INDEPENDENT_CODE ON)

add_library(phasefields::cost ALIAS cost)

