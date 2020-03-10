cmake_minimum_required(VERSION 3.12)

project(decrypt_core)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost COMPONENTS program_options REQUIRED)

find_package(SQLite3 REQUIRED)

find_package(GDAL REQUIRED)

include_directories(headeronly_src include)

add_executable(model_1 model_1.cpp)
add_executable(model_2 model_2.cpp)

target_link_libraries(model_1 LINK_PUBLIC ${SQLite3_LIBRARIES} GDAL::GDAL ${Boost_LIBRARIES})
target_link_libraries(model_2 LINK_PUBLIC ${SQLite3_LIBRARIES} GDAL::GDAL ${Boost_LIBRARIES})

# INSTALL(TARGETS spatial_process1 DESTINATION ${STAGED_INSTALL_PREFIX}/${PROJECT_NAME})
install(
  TARGETS
    model_1 model_2
  RUNTIME
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT bin
  )
