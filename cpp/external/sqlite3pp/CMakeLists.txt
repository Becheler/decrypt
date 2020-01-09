cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(
  "sqlite3pp"
  LANGUAGES CXX
  DESCRIPTION "SQLite3 API more friendly to C++ users"
  HOMEPAGE_URL "https://github.com/iwongu/sqlite3pp")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(${PROJECT_NAME} INTERFACE)
# add alias so the project can be used with add_subdirectory
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

target_include_directories(${PROJECT_NAME} INTERFACE headeronly_src)
