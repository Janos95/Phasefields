set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

add_library(Mesh STATIC
        Mesh.cpp
        Mesh.h
        MeshElements.cpp
        MeshElements.h
        MeshFeature.cpp
        MeshFeature.h
        FastMarchingMethod.cpp
        FastMarchingMethod.h
        Surface.h
        Dijkstra.h
        GraphCommon.h
        FEM.cpp
        FEM.h
        )

target_include_directories(Mesh PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_link_libraries(Mesh PRIVATE
        ScopedTimer::ScopedTimer
        Corrade::Containers
        Magnum::GL
        #Phasefield::Optimization
        )

target_link_libraries(Mesh PUBLIC
        Phasefield::Utilities
        )

set_property(TARGET Mesh PROPERTY POSITION_INDEPENDENT_CODE ON)

add_library(Phasefield::Mesh ALIAS Mesh)
