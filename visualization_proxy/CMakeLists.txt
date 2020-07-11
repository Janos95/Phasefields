add_library(visualization_proxy STATIC
        visualization_proxy.cpp
        visualization_proxy.hpp
        )

target_include_directories(visualization_proxy PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(visualization_proxy
        PRIVATE
            Magnum::GL
            Magnum::Magnum
            phasefields::viewer
        PUBLIC
            phasefields::utilities
        )

set_property(TARGET visualization_proxy PROPERTY POSITION_INDEPENDENT_CODE ON)
add_library(phasefields::visualization_proxy ALIAS visualization)
