add_library(utilities INTERFACE)
target_sources(utilities INTERFACE
        unique_function.h
        detail/pointer_union.h
        detail/pointer_int_pair.h
        detail/pointer_like_type_traits.h
        shared_ressource.hpp
        function_ref.hpp
        normalizeInto.cpp
        normalizeInto.hpp
        smart_enum.hpp
        permutations.hpp
        stl_algorithm.hpp
        y_combinator.hpp
        types.hpp
        )

target_include_directories(utilities INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
add_library(phasefields::utilities ALIAS utilities)

#add_executable(bench hash_map_bench.cpp)
#target_link_libraries(bench PUBLIC utilities Magnum::Magnum Eigen3::Eigen)

#add_executable(map_test test.cpp)