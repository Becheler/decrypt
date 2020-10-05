cmake_minimum_required(VERSION 3.12)

project(decrypt_core)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost COMPONENTS program_options REQUIRED)

find_package(SQLite3 REQUIRED)
IF (SQLite3_FOUND)
    INCLUDE_DIRECTORIES(${SQLite3_INCLUDE_DIRS})
ENDIF()

find_package(GDAL REQUIRED)
IF (GDAL_FOUND)
    INCLUDE_DIRECTORIES(${GDAL_INCLUDE_DIR})
ENDIF()

include_directories(headeronly_src include)

add_executable(model_1 model_1.cpp)
add_executable(model_2 model_2.cpp)
add_executable(model_3 model_3.cpp)

target_link_libraries(model_1 LINK_PUBLIC ${SQLite3_LIBRARIES} GDAL::GDAL ${Boost_LIBRARIES})
target_link_libraries(model_2 LINK_PUBLIC ${SQLite3_LIBRARIES} GDAL::GDAL ${Boost_LIBRARIES})
target_link_libraries(model_3 LINK_PUBLIC ${SQLite3_LIBRARIES} GDAL::GDAL ${Boost_LIBRARIES})

add_executable(model_1_data model_1_data.cpp)
target_link_libraries(model_1_data LINK_PUBLIC ${SQLite3_LIBRARIES} GDAL::GDAL ${Boost_LIBRARIES})

add_executable(model_3_data model_3_data.cpp)
target_link_libraries(model_3_data LINK_PUBLIC ${SQLite3_LIBRARIES} GDAL::GDAL ${Boost_LIBRARIES})

# INSTALL(TARGETS spatial_process1 DESTINATION ${STAGED_INSTALL_PREFIX}/${PROJECT_NAME})
install(
  TARGETS
    model_1 model_2 model_3 model_1_data model_3_data
  RUNTIME
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT bin
  )
