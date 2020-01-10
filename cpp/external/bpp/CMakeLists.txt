cmake_minimum_required(VERSION 3.12)

project("bpp" LANGUAGES C)

add_executable(${PROJECT_NAME} src/bpp.c src/tools.c)
target_include_directories(${PROJECT_NAME} INTERFACE src)
target_link_libraries(${PROJECT_NAME} m)
target_compile_options(${PROJECT_NAME} PRIVATE -O3)
set_property(TARGET ${PROJECT_NAME} PROPERTY C_STANDARD 99)
