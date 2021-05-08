#set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

#set(CUDA_TOOLKIT_ROOT_DIR /opt/cuda)
#set(CMAKE_CUDA_COMPILER /opt/cuda/bin/nvcc)
#enable_language(CUDA)
#
#add_library(CudaCG
#        CudaCG.cpp
#        CudaCG.h
#        )
#
#find_package(OpenMP)
#find_package(CUDAToolkit)
#
#target_link_libraries(CudaCG PUBLIC CUDA::cusparse CUDA::cusolver CUDA::cudart Eigen3::Eigen ScopedTimer::ScopedTimer)
#add_library(Phasefield::CudaCG ALIAS CudaCG)

#add_library(SparseMatrix
#        SparseMatrix.cpp
#        SparseMatrix.h
#        )
#
#add_library(Phasefield::SparseMatrix ALIAS SparseMatrix)

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
        LinearizedElasticity.cpp
        LinearizedElasticity.h
        )

target_include_directories(Cost PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_include_directories(Cost SYSTEM PRIVATE ${Eigen3_INCLUDE_DIRS})

target_link_libraries(Cost PRIVATE
        ScopedTimer::ScopedTimer
        Phasefield::Viewer
        #Phasefield::CudaCG
        )

find_package(OpenMP)

target_link_libraries(Cost PUBLIC
        Phasefield::Utilities
        #Phasefield::Visualization
        Phasefield::VisualizationProxy
        Phasefield::Mesh
        Phasefield::Optimization
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



