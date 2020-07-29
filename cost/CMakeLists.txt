find_package(Ceres CONFIG REQUIRED)
find_package(TBB CONFIG REQUIRED)

find_package(PkgConfig REQUIRED)

list(APPEND CMAKE_PREFIX_PATH /home/janos/libraries/prefix/)
set(PKG_CONFIG_USE_CMAKE_PREFIX_PATH ON)

pkg_check_modules(Adolc REQUIRED IMPORTED_TARGET adolc)

set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

add_library(cost STATIC
        fem.cpp
        fem.hpp
        modica_mortola.cpp
        modica_mortola.hpp
        strang_rules.hpp)

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


