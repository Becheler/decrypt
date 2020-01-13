cmake_minimum_required(VERSION 3.12)

project(decrypt)

enable_testing()

add_subdirectory(cpp)

INSTALL(DIRECTORY assets/ python/ R/ example DESTINATION ${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME})
