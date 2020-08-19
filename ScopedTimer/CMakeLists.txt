
add_library(ScopedTimer STATIC include/ScopedTimer/ScopedTimer.cpp)

target_include_directories(ScopedTimer PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

target_link_libraries(ScopedTimer PUBLIC Corrade::Containers)
set_property(TARGET ScopedTimer PROPERTY POSITION_INDEPENDENT_CODE ON)
add_library(ScopedTimer::ScopedTimer ALIAS ScopedTimer)