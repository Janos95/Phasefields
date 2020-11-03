add_library(Utilities
        UniqueFunction.h
        SharedPointer.h
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
        Utility.h
        Traits.h
        Algorithms.h
        )

target_include_directories(Utilities PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(Utilities PUBLIC Corrade::Containers Magnum::Magnum)
add_library(Phasefield::Utilities ALIAS Utilities)

#add_executable(bench hash_map_bench.cpp)
#target_link_libraries(bench PUBLIC utilities Magnum::Magnum Eigen3::Eigen)
