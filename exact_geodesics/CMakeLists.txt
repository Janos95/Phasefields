

add_library(exact_geodesics INTERFACE)

target_sources(exact_geodesics INTERFACE
        geodesic_algorithm_graph_base.h
        geodesic_algorithm_subdivision.h
        geodesic_constants_and_simple_functions.h
        geodesic_memory.h
        geodesic_mesh.h
        geodesic_algorithm_base.h
        geodesic_algorithm_dijkstra.h
        geodesic_algorithm_dijkstra_alternative.h
        geodesic_algorithm_exact.h
        geodesic_algorithm_exact_elements.h
        geodesic_mesh_elements.h
        )

target_include_directories(exact_geodesics SYSTEM INTERFACE ${CMAKE_CURRENT_LIST_DIR})

add_library(phasefields::exact_geodesics ALIAS exact_geodesics)
