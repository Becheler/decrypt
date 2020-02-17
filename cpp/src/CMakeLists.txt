cmake_minimum_required(VERSION 3.12)

project(decrypt)

set(CMAKE_CXX_STANDARD 17)

FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)

FIND_PACKAGE(GDAL)
IF (GDAL_FOUND)
    INCLUDE_DIRECTORIES(${GDAL_INCLUDE_DIR})
ENDIF()

find_package (SQLite3 REQUIRED)

include_directories(headeronly_src include)

add_executable(spatial_process1 main.cpp)
target_link_libraries(spatial_process1 LINK_PUBLIC ${SQLite3_LIBRARIES} ${GDAL_LIBRARY} ${Boost_LIBRARIES})
#INSTALL(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_SOURCE_DIR}/assets)
#INSTALL(TARGETS
#INSTALL(TARGETS # Creates a folder called "executables" and adds target
# project (helloWorld.vcproj) under it

INSTALL(TARGETS spatial_process1 DESTINATION ${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME})
