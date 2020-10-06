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
        DiffuseYamabe.cpp
        DiffuseYamabe.h
        )

target_include_directories(Cost PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_link_libraries(Cost PRIVATE
        ScopedTimer::ScopedTimer
        Eigen3::Eigen
        Phasefield::Viewer
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


if (NOT CORRADE_TARGET_EMSCRIPTEN)
    target_link_libraries(Cost PRIVATE
            OpenMP::OpenMP_CXX
            TBB::tbb
            PkgConfig::Adolc
            SuiteSparse::UMFPACK
            )

endif ()


