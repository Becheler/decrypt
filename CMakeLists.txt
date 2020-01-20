cmake_minimum_required(VERSION 3.12)

project(decrypt)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

message(STATUS: "Buildtype: ${CMAKE_BUILD_TYPE}")

message(STATUS "C++ flags, Debug configuration: ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "C++ flags, Release configuration: ${CMAKE_CXX_FLAGS_RELEASE}")

enable_testing()

add_subdirectory(cpp)

INSTALL(DIRECTORY assets/ python/ R/ example DESTINATION ${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME})
