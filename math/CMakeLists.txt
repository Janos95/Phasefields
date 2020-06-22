list(APPEND CMAKE_PREFIX_PATH ${CMAKE_CURRENT_LIST_DIR}/openblas_module)
find_package(OpenBLAS CONFIG)

add_library(math STATIC
        detail.hpp
        tensor.cpp
        tensor.hpp
        matrix.hpp
        matrix_expressions.hpp
        tensor_expressions.hpp
        decomposition.hpp
        )

set_property(TARGET math PROPERTY CXX_STANDARD 20)
target_link_libraries(math PUBLIC OpenBLAS::OpenBLAS)