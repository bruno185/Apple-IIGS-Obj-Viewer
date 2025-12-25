cmake_minimum_required(VERSION 3.10)
project(GS3DpViewer C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

find_package(SDL2 CONFIG)
if (SDL2_FOUND)
    set(SDL2_LIBS ${SDL2_LIBRARIES})
else()
    if (NOT MSVC)
        find_package(PkgConfig REQUIRED)
        set(_pkg_ok FALSE)
        if (PKG_CONFIG_FOUND)
            pkg_check_modules(SDL2_PKG QUIET sdl2)
            if (SDL2_PKG_FOUND)
                set(_pkg_ok TRUE)
                set(SDL2_INCLUDE_DIRS ${SDL2_PKG_INCLUDE_DIRS})
                set(SDL2_LIBS ${SDL2_PKG_LIBRARIES})
            endif()
        endif()
        if (NOT _pkg_ok)
            # Fallback: try common MinGW locations
            find_path(SDL2_INCLUDE_DIR SDL.h PATHS /mingw64/include /usr/include /usr/local/include C:/msys64/mingw64/include)
            find_library(SDL2_LIBRARY NAMES SDL2 SDL2main PATHS /mingw64/lib /usr/lib /usr/local/lib C:/msys64/mingw64/lib)
            if (SDL2_INCLUDE_DIR AND SDL2_LIBRARY)
                set(SDL2_INCLUDE_DIRS ${SDL2_INCLUDE_DIR})
                set(SDL2_LIBS ${SDL2_LIBRARY})
            endif()
        endif()
    endif()
endif()

# Console SDL viewer (skip on native Windows build to avoid confusion)
if(NOT WIN32)
add_executable(viewer
    src/main.c
    src/fixed32.h
    src/obj.c
    src/painter.c
    src/renderer.c
)

target_include_directories(viewer PRIVATE src ${SDL2_INCLUDE_DIRS})

target_link_libraries(viewer ${SDL2_LIBS})

target_compile_definitions(viewer PRIVATE USE_SDL)
endif()
if (WIN32)
    add_executable(viewer_win32 WIN32
        src/win32_main.c
        src/win32_renderer.c
        src/fixed32.h
        src/obj.c
        src/painter.c
        src/viewer.h
    )
    target_include_directories(viewer_win32 PRIVATE src)
    target_link_libraries(viewer_win32 gdi32 user32)
endif()
