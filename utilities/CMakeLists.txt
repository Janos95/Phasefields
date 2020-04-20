
add_library(utilities INTERFACE)
target_sources(utilities INTERFACE
        robin_map.h
        robin_growth_policy.h
        robin_hash.h
        unique_function.h
        pointer_union.h
        pointer_int_pair.h
        pointer_like_type_traits.h)
target_include_directories(utilities INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
add_library(phasefields::utilities ALIAS utilities)