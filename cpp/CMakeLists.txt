cmake_minimum_required(VERSION 3.12)

project(decrypt)

set(CMAKE_CXX_STANDARD 17)
add_compile_options(-O3 -Wall -pedantic-errors -Wextra -DNDEBUG)

find_package(QUETZAL CONFIG REQUIRED)

FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
ENDIF()


FIND_PACKAGE(GDAL)
IF (GDAL_FOUND)
    INCLUDE_DIRECTORIES(${GDAL_INCLUDE_DIR})
ENDIF()


find_package(SQLite3)
include_directories(${SQLite3_INCLUDE_DIRS})


find_path(SQLITE3PP_INCLUDE_DIR sqlite3pp.h)
if((NOT SQLITE3PP_INCLUDE_DIR) OR (NOT EXISTS ${SQLITE3PP_INCLUDE_DIR}))
  message("Unable to find sqlite3pp")
  execute_process(COMMAND git submodule update --init -- external/sqlite3pp
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  set(SQLITE3PP_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/sqlite3pp/headeronly_src CACHE PATH "sqlite3pp include directory")
  # also install it
  install(DIRECTORY ${SQLITE3PP_INCLUDE_DIR} DESTINATION /home/becheler/dev/libs/sqlite3pp)
  # for convenience setup a target
  add_library(sqlite3pp INTERFACE)
  target_include_directories(sqlite3pp INTERFACE
                             $<BUILD_INTERFACE:${SQLITE3PP_INCLUDE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

  target_include_directories(
   sqlite3pp
   INTERFACE $<BUILD_INTERFACE:${sqlite3pp_SOURCE_DIR}/include>
             $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
  # need to export target as well
  install(TARGETS sqlite3pp EXPORT my_export_set DESTINATION /home/becheler/dev/libs/sqlite3pp)
else()
  message("found sqlite3pp")
  add_library(sqlite3pp INTERFACE)
	target_include_directories(sqlite3pp INTERFACE ${SQLITE3PP_INCLUDE_DIR})
	target_link_libraries(sqlite3pp INTERFACE ${SQLITE3PP_LIBRARY})
endif()

#include_directories( /home/becheler/dev/)

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} LINK_PUBLIC QUETZAL::QUETZAL ${GDAL_LIBRARY}  ${SQLITE3_LIBRARY} ${Boost_LIBRARIES} ${SQLITE3PP_LIBRARY} )
