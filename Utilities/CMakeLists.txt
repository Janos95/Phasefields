add_library(Utilities
        UniqueFunction.h
        detail/pointer_union.h
        detail/pointer_int_pair.h
        detail/pointer_like_type_traits.h
        SharedRessource.h
        FunctionRef.h
        normalizeInto.cpp
        normalizeInto.hpp
        SmartEnum.h
        StlAlgorithm.h
        YCombinator.h
        c1_functions.hpp
        Enums.h
        Heap.h
        CircularBuffer.h
        )

target_include_directories(Utilities PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(Utilities PUBLIC Corrade::Containers)
add_library(Phasefield::Utilities ALIAS Utilities)

#add_executable(bench hash_map_bench.cpp)
#target_link_libraries(bench PUBLIC utilities Magnum::Magnum Eigen3::Eigen)
