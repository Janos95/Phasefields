find_package(Ceres CONFIG REQUIRED)
find_package(TBB CONFIG REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(Ipopt REQUIRED IMPORTED_TARGET ipopt)
pkg_check_modules(Adolc REQUIRED IMPORTED_TARGET adolc)

file(GLOB cost_headers ${CMAKE_CURRENT_SOURCE_DIR} *.hpp)
file(GLOB cost_sources ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

add_library(cost STATIC ${cost_sources} ${cost_headers})
target_include_directories(cost PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

find_package(libigl CONFIG REQUIRED)

target_link_libraries(cost PRIVATE
        ceres
        igl::core
        Eigen3::Eigen
        autodiff
        TBB::tbb
        phasefields::scoped_timer
        )

target_link_libraries(cost PUBLIC
        PkgConfig::Ipopt
        PkgConfig::Adolc
        phasefields::utilities
        )

add_library(phasefields::cost ALIAS cost)

