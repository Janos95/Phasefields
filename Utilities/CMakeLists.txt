add_library(Utilities
        UniqueFunction.h
        detail/pointer_union.h
        detail/pointer_int_pair.h
        detail/pointer_like_type_traits.h
        SharedRessource.h
        FunctionRef.h
        SmartEnum.h
        StlAlgorithm.h
        YCombinator.h
        C1Functions.h
        Enums.h
        Heap.h
        CircularBuffer.h
        SparseMatrix.cpp
        SparseMatrix.h
        )

target_include_directories(Utilities PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(Utilities PUBLIC Corrade::Containers Magnum::Magnum)
add_library(Phasefield::Utilities ALIAS Utilities)

#add_executable(bench hash_map_bench.cpp)
#target_link_libraries(bench PUBLIC utilities Magnum::Magnum Eigen3::Eigen)
