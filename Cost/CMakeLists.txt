find_package(Ceres CONFIG REQUIRED)
find_package(TBB CONFIG REQUIRED)

find_package(PkgConfig REQUIRED)

list(APPEND CMAKE_PREFIX_PATH /home/janos/libraries/prefix/)
set(PKG_CONFIG_USE_CMAKE_PREFIX_PATH ON)

pkg_check_modules(Adolc REQUIRED IMPORTED_TARGET adolc)

set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

add_library(Cost STATIC

        #Fem.cpp
        #Fem.h
        #ModicaMortola.cpp
        #ModicaMortola.h
        #SparseMatrix.h
        #SparseMatrix.cpp
        #StrangRules.hpp
        #FastMarchingMethod.cpp
        #FastMarchingMethod.h
        #YamabeSolver.cpp
        #YamabeSolver.h
        )

target_include_directories(Cost PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_link_libraries(Cost PRIVATE
        ScopedTimer::ScopedTimer
        #Phasefield::Optimization
        )

target_link_libraries(Cost PUBLIC
        Phasefield::Utilities
        Phasefield::Visualization
        Phasefield::VisualizationProxy
        Eigen3::Eigen
        )

set_property(TARGET Cost PROPERTY POSITION_INDEPENDENT_CODE ON)

add_library(Phasefield::Cost ALIAS Cost)


