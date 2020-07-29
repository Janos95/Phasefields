
set(OpenBLAS_INCLUDE_DIRS_RELEASE "/home/janos/.conan/data/openblas/0.3.9/_/_/package/141e6f5bf5b83e010f9dc6c2cce61203915d16ec/include")
set(OpenBLAS_INCLUDE_DIR_RELEASE "/home/janos/.conan/data/openblas/0.3.9/_/_/package/141e6f5bf5b83e010f9dc6c2cce61203915d16ec/include")
set(OpenBLAS_INCLUDES_RELEASE "/home/janos/.conan/data/openblas/0.3.9/_/_/package/141e6f5bf5b83e010f9dc6c2cce61203915d16ec/include")
set(OpenBLAS_RES_DIRS_RELEASE)
set(OpenBLAS_DEFINITIONS_RELEASE)
set(OpenBLAS_LINKER_FLAGS_RELEASE_LIST "" "")
set(OpenBLAS_COMPILE_DEFINITIONS_RELEASE)
set(OpenBLAS_COMPILE_OPTIONS_RELEASE_LIST "" "")
set(OpenBLAS_LIBRARIES_TARGETS_RELEASE "") # Will be filled later, if CMake 3
set(OpenBLAS_LIBRARIES_RELEASE "") # Will be filled later
set(OpenBLAS_LIBS_RELEASE "") # Same as OpenBLAS_LIBRARIES
set(OpenBLAS_SYSTEM_LIBS_RELEASE pthread)
set(OpenBLAS_FRAMEWORK_DIRS_RELEASE)
set(OpenBLAS_FRAMEWORKS_RELEASE)
set(OpenBLAS_FRAMEWORKS_FOUND_RELEASE "") # Will be filled later
set(OpenBLAS_BUILD_MODULES_PATHS_RELEASE)

conan_find_apple_frameworks(OpenBLAS_FRAMEWORKS_FOUND_RELEASE "${OpenBLAS_FRAMEWORKS_RELEASE}" "${OpenBLAS_FRAMEWORK_DIRS_RELEASE}")

mark_as_advanced(OpenBLAS_INCLUDE_DIRS_RELEASE
        OpenBLAS_INCLUDE_DIR_RELEASE
        OpenBLAS_INCLUDES_RELEASE
        OpenBLAS_DEFINITIONS_RELEASE
        OpenBLAS_LINKER_FLAGS_RELEASE_LIST
        OpenBLAS_COMPILE_DEFINITIONS_RELEASE
        OpenBLAS_COMPILE_OPTIONS_RELEASE_LIST
        OpenBLAS_LIBRARIES_RELEASE
        OpenBLAS_LIBS_RELEASE
        OpenBLAS_LIBRARIES_TARGETS_RELEASE)

# Find the real .lib/.a and add them to OpenBLAS_LIBS and OpenBLAS_LIBRARY_LIST
set(OpenBLAS_LIBRARY_LIST_RELEASE openblas)
set(OpenBLAS_LIB_DIRS_RELEASE "/home/janos/.conan/data/openblas/0.3.9/_/_/package/141e6f5bf5b83e010f9dc6c2cce61203915d16ec/lib")

# Gather all the libraries that should be linked to the targets (do not touch existing variables):
set(_OpenBLAS_DEPENDENCIES_RELEASE "${OpenBLAS_FRAMEWORKS_FOUND_RELEASE} ${OpenBLAS_SYSTEM_LIBS_RELEASE} ")

conan_package_library_targets("${OpenBLAS_LIBRARY_LIST_RELEASE}"  # libraries
        "${OpenBLAS_LIB_DIRS_RELEASE}"      # package_libdir
        "${_OpenBLAS_DEPENDENCIES_RELEASE}"  # deps
        OpenBLAS_LIBRARIES_RELEASE            # out_libraries
        OpenBLAS_LIBRARIES_TARGETS_RELEASE    # out_libraries_targets
        "_RELEASE"                          # build_type
        "OpenBLAS")                                      # package_name

set(OpenBLAS_LIBS_RELEASE ${OpenBLAS_LIBRARIES_RELEASE})

foreach (_FRAMEWORK ${OpenBLAS_FRAMEWORKS_FOUND_RELEASE})
    list(APPEND OpenBLAS_LIBRARIES_TARGETS_RELEASE ${_FRAMEWORK})
    list(APPEND OpenBLAS_LIBRARIES_RELEASE ${_FRAMEWORK})
endforeach ()

foreach (_SYSTEM_LIB ${OpenBLAS_SYSTEM_LIBS_RELEASE})
    list(APPEND OpenBLAS_LIBRARIES_TARGETS_RELEASE ${_SYSTEM_LIB})
    list(APPEND OpenBLAS_LIBRARIES_RELEASE ${_SYSTEM_LIB})
endforeach ()

# We need to add our requirements too
set(OpenBLAS_LIBRARIES_TARGETS_RELEASE "${OpenBLAS_LIBRARIES_TARGETS_RELEASE};")
set(OpenBLAS_LIBRARIES_RELEASE "${OpenBLAS_LIBRARIES_RELEASE};")

set(CMAKE_MODULE_PATH "/home/janos/.conan/data/openblas/0.3.9/_/_/package/141e6f5bf5b83e010f9dc6c2cce61203915d16ec/" ${CMAKE_MODULE_PATH})
set(CMAKE_PREFIX_PATH "/home/janos/.conan/data/openblas/0.3.9/_/_/package/141e6f5bf5b83e010f9dc6c2cce61203915d16ec/" ${CMAKE_PREFIX_PATH})

foreach (_BUILD_MODULE_PATH ${OpenBLAS_BUILD_MODULES_PATHS_RELEASE})
    include(${_BUILD_MODULE_PATH})
endforeach ()
