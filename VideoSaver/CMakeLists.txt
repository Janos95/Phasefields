
add_library(VideoSaver STATIC VideoSaver.cpp VideoSaver.h)

target_include_directories(VideoSaver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
add_library(Phasefield::VideoSaver ALIAS VideoSaver)

if (NOT CORRADE_TARGET_EMSCRIPTEN)
    target_link_libraries(VideoSaver
            PUBLIC
            PkgConfig::libavcodec
            PkgConfig::libavutil
            PkgConfig::libswscale
            Corrade::Containers
            Phasefield::Utilities
            )
endif ()

