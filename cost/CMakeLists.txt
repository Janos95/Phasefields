find_package(Ceres CONFIG REQUIRED)
find_package(TBB CONFIG REQUIRED)

find_package(PkgConfig REQUIRED)

list(APPEND CMAKE_PREFIX_PATH /home/janos/libraries/prefix/)
set(PKG_CONFIG_USE_CMAKE_PREFIX_PATH ON)

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
        graph_common.hpp
        meta_data.hpp
        meta_data.cpp
        modica_mortola.cpp
        modica_mortola.hpp
        stopping_criteria.cpp
        stopping_criteria.hpp
        union_find.cpp
        union_find.hpp
        strang_rules.hpp
        phasefield_tree.hpp
        )

target_include_directories(cost PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_link_libraries(cost PRIVATE
        phasefields::scoped_timer
        phasefields::optimization
        )

target_link_libraries(cost PUBLIC
        Eigen3::Eigen
        phasefields::utilities
        phasefields::visualization
        phasefields::visualization_proxy
        )

set_property(TARGET cost PROPERTY POSITION_INDEPENDENT_CODE ON)

add_library(phasefields::cost ALIAS cost)


