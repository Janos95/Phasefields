find_package(PkgConfig REQUIRED)
pkg_check_modules(libavcodec REQUIRED IMPORTED_TARGET libavcodec)
pkg_check_modules(libavutil REQUIRED IMPORTED_TARGET libavutil)
pkg_check_modules(libswscale REQUIRED IMPORTED_TARGET libswscale)

add_library(VideoSaver STATIC VideoSaver.cpp VideoSaver.h)
target_link_libraries(VideoSaver
        PUBLIC
            PkgConfig::libavcodec
            PkgConfig::libavutil
            PkgConfig::libswscale
            Corrade::Containers
            Phasefield::Utilities
        )

target_include_directories(VideoSaver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
add_library(Phasefield::VideoSaver ALIAS VideoSaver)
