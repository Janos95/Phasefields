set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

set(CUDA_TOOLKIT_ROOT_DIR /opt/cuda)
set(CMAKE_CUDA_COMPILER /opt/cuda/bin/nvcc)
enable_language(CUDA)
add_library(CudaSolver
        CudaSolver.cpp
        CudaSolver.h
        CuSparse/cusparse_cholesky_solver.cpp
        CuSparse/cusparse_cholesky_solver.h
        CuSparse/cusparse_wrapper.cpp
        CuSparse/cusparse_wrapper.h
        )

find_package(OpenMP)
find_package(CUDAToolkit)

target_link_libraries(CudaSolver PUBLIC Eigen3::Eigen CUDA::cusparse CUDA::cusolver CUDA::cudart)
add_library(Phasefield::CudaSolver ALIAS CudaSolver)


add_library(Cost STATIC
        Cost.h
        ModicaMortola.cpp
        ModicaMortola.h
        ConnectednessConstraint.h
        ConnectednessConstraint.cpp
        #StrangRules.hpp
        CudaSolver.cpp
        CudaSolver.h
        StoppingCriteria.cpp
        StoppingCriteria.h
        UnionFind.cpp
        UnionFind.h
        DiffuseYamabe.cpp
        DiffuseYamabe.h
        LinearizedElasticity.cpp
        LinearizedElasticity.h
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
        Phasefield::CudaSolver
        OpenMP::OpenMP_CXX
        )

set_property(TARGET Cost PROPERTY POSITION_INDEPENDENT_CODE ON)
add_library(Phasefield::Cost ALIAS Cost)

if (PHASEFIELD_WITH_ADOLC)
    target_link_libraries(Cost PRIVATE PkgConfig::Adolc)
    target_compile_definitions(Cost PRIVATE PHASEFIELD_WITH_ADOLC)
endif()

if (PHASEFIELD_WITH_TBB)
    target_compile_definitions(Cost PRIVATE PHASEFIELD_WITH_TBB)
endif()

if (PHASEFIELD_WITH_SUITESPARSE)
    target_link_libraries(Cost PRIVATE SuiteSparse::UMFPACK)
    target_compile_definitions(Cost PRIVATE PHASEFIELD_WITH_SUITESPARSE)
endif()



