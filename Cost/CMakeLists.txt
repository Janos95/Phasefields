find_package(Ceres CONFIG REQUIRED)
find_package(TBB CONFIG REQUIRED)
find_package(OpenMP)
find_package(PkgConfig REQUIRED)

list(APPEND CMAKE_PREFIX_PATH /home/janos/libraries/prefix/)
set(PKG_CONFIG_USE_CMAKE_PREFIX_PATH ON)

pkg_check_modules(Adolc REQUIRED IMPORTED_TARGET adolc)

set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

add_library(Cost STATIC
        Cost.h
        ModicaMortola.cpp
        ModicaMortola.h
        ConnectednessConstraint.h
        ConnectednessConstraint.cpp
        #StrangRules.hpp
        StoppingCriteria.cpp
        StoppingCriteria.h
        UnionFind.cpp
        UnionFind.h
        WeakYamabe.cpp
        WeakYamabe.h
        )

target_include_directories(Cost PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_link_libraries(Cost PRIVATE
        ScopedTimer::ScopedTimer
        OpenMP::OpenMP_CXX
        Eigen3::Eigen
        ImGui::ImGui
        Phasefield::Viewer
        TBB::tbb
        PkgConfig::Adolc
        )

target_link_libraries(Cost PUBLIC
        Phasefield::Utilities
        #Phasefield::Visualization
        Phasefield::VisualizationProxy
        Phasefield::Mesh
        Phasefield::Optimization
        )

set_property(TARGET Cost PROPERTY POSITION_INDEPENDENT_CODE ON)

add_library(Phasefield::Cost ALIAS Cost)


