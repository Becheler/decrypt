cmake_minimum_required(VERSION 3.12)

project(spatial_process)

set(CMAKE_CXX_STANDARD 17)
#add_compile_options(-O3 -Wall -pedantic-errors -Wextra -DNDEBUG)

FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
ENDIF()

FIND_PACKAGE(GDAL)
IF (GDAL_FOUND)
    INCLUDE_DIRECTORIES(${GDAL_INCLUDE_DIR})
ENDIF()

FIND_PACKAGE(SQLite3)
IF (SQLite3_FOUND)
    INCLUDE_DIRECTORIES(${SQLite3_INCLUDE_DIRS})
ENDIF()

enable_testing()

add_subdirectory(external)

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} LINK_PUBLIC quetzal sqlite3pp ${GDAL_LIBRARY}  ${SQLite3_LIBRARIES} ${Boost_LIBRARIES})
INSTALL(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_SOURCE_DIR}/assets)
