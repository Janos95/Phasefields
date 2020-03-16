file(GLOB cost_headers ${CMAKE_CURRENT_SOURCE_DIR} *.hpp)
file(GLOB cost_sources ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

add_library(cost STATIC ${cost_sources} ${cost_headers})
target_include_directories(cost PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_link_libraries(cost PUBLIC ceres)
target_link_libraries(cost PRIVATE igl::core Folly::folly autodiff phasefields::scoped_timer)
if(OpenMP_CXX_FOUND)
    target_link_libraries(cost PRIVATE OpenMP::OpenMP_CXX)
endif()
add_library(phasefields::cost ALIAS cost)

