cmake_minimum_required(VERSION 3.13)
file(GLOB_RECURSE SOURCES RELATIVE ${CMAKE_SOURCE_DIR} "libft/*.c")
file(GLOB_RECURSE HEADERS RELATIVE ${CMAKE_SOURCE_DIR} "libft/includes/*.h")
file(GLOB_RECURSE SOURCES1 RELATIVE ${CMAKE_SOURCE_DIR} "libgm/*.c")
file(GLOB_RECURSE HEADERS1 RELATIVE ${CMAKE_SOURCE_DIR} "libgm/include/*.h")
file(GLOB_RECURSE SOURCES2 RELATIVE ${CMAKE_SOURCE_DIR} "src/*.c")
file(GLOB_RECURSE HEADERS2 RELATIVE ${CMAKE_SOURCE_DIR} "includes/*.h")
project(RTv1 C)

set(CMAKE_C_STANDARD 99)

#add_compile_options(-Wall -Werror -Wextra)

#target_link_libraries(RTv1 "-pthread")
include_directories("includes")
include_directories("libft/includes")
link_directories(/usr/local/lib)
link_libraries(mlx)
find_library(OPEN_GL OpenGL)
find_library(APP_KIT AppKit)

add_executable(RTv1 ${SOURCES} ${HEADERS} ${HEADERS1} ${SOURCES1} ${HEADERS2} ${SOURCES2} src/cone.c src/cylinder.c src/plane.c)
target_link_libraries(RTv1 ${OPEN_GL} ${APP_KIT})

cmake_minimum_required(VERSION 3.14)
