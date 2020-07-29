
corrade_add_resource(TMO_Rcs resources.conf)

add_library(shaders SHARED
        Phong.h
        Phong.cpp
        TextureCoordinates.h
        TextureCoordinates.cpp
        ${TMO_Rcs} TextureCoordinates.cpp TextureCoordinates.h ScreenTriangle.cpp ScreenTriangle.h Remap.cpp Remap.h Combine.cpp Combine.h Sobel.cpp Sobel.h Reduction.cpp Reduction.h ChainRule.cpp ChainRule.h)


target_include_directories(shaders PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(shaders
        PUBLIC
        Magnum::Magnum
        Magnum::Shaders
        scoped_timer::scoped_timer
        )

set_property(TARGET shaders PROPERTY POSITION_INDEPENDENT_CODE ON)
add_library(tmo::shaders ALIAS shaders)