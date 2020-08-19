
corrade_add_resource(Phasefields_Rcs resources.conf)

add_library(Shaders SHARED
        Phong.h
        Phong.cpp
        ${Phasefields_Rcs})


target_include_directories(Shaders PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(Shaders
        PUBLIC
        Magnum::Magnum
        Magnum::Shaders
        ScopedTimer::ScopedTimer
        )

set_property(TARGET Shaders PROPERTY POSITION_INDEPENDENT_CODE ON)
add_library(Phasefield::Shaders ALIAS Shaders)