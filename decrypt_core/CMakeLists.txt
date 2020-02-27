cmake_minimum_required(VERSION 3.12)

project(decrypt_core)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost COMPONENTS program_options REQUIRED)

find_package(SQLite3 REQUIRED)

find_package(GDAL REQUIRED)

include_directories(headeronly_src include)

add_executable(spatial_process1 main.cpp)

target_link_libraries(spatial_process1 LINK_PUBLIC ${SQLite3_LIBRARIES} GDAL::GDAL ${Boost_LIBRARIES})

# INSTALL(TARGETS spatial_process1 DESTINATION ${STAGED_INSTALL_PREFIX}/${PROJECT_NAME})
install(
  TARGETS
    spatial_process1
  RUNTIME
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT bin
  )
