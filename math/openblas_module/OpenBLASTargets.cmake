
if(NOT TARGET OpenBLAS::OpenBLAS)
    add_library(OpenBLAS::OpenBLAS INTERFACE IMPORTED)
endif()

# Load the debug and release library finders
get_filename_component(_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
file(GLOB CONFIG_FILES "${_DIR}/OpenBLASTarget-*.cmake")

foreach(f ${CONFIG_FILES})
  include(${f})
endforeach()
    
