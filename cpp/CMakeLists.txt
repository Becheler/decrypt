cmake_minimum_required(VERSION 2.8)

project(main)

set(CMAKE_CXX_STANDARD 17)
add_compile_options(-O3 -Wall -pedantic-errors -Wextra -DNDEBUG -lsqlite3)

FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
ENDIF()

FIND_PACKAGE(GDAL)
IF (GDAL_FOUND)
    INCLUDE_DIRECTORIES(${GDAL_INCLUDE_DIR})
ENDIF()

find_path(SQLITE3_INCLUDE_DIR sqlite3.h)
find_library(SQLITE3_LIBRARY sqlite3)
if(NOT SQLITE3_INCLUDE_DIR)
  message(SEND_ERROR "Failed to find sqlite3")
endif()
if(NOT SQLITE3_LIBRARY)
  message(SEND_ERROR "Failed to find the sqlite3 library")
endif()

include_directories( /home/becheler/dev/)

add_executable(main main.cpp)
TARGET_LINK_LIBRARIES(main LINK_PUBLIC ${GDAL_LIBRARY}  ${SQLITE3_LIBRARY} ${Boost_LIBRARIES} )
