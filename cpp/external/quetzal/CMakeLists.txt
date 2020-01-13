cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(
  "quetzal"
  VERSION 0.1
  LANGUAGES CXX
  DESCRIPTION "A header only c++ library for spatially explicit coalescence simulation."
  HOMEPAGE_URL "https://github.com/becheler/quetzal")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message(STATUS "Project will be installed to ${CMAKE_INSTALL_PREFIX}")

include(GNUInstallDirs)

add_library(${PROJECT_NAME} INTERFACE)
# add alias so the project can be used with add_subdirectory
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

# Include Boost as an imported target
find_package(Boost 1.61 REQUIRED COMPONENTS unit_test_framework filesystem)

# Include Gdal
find_package(GDAL REQUIRED)
IF (GDAL_FOUND)
    INCLUDE_DIRECTORIES(${GDAL_INCLUDE_DIR})
ENDIF()

# Adding the install interface generator expression makes sure that the include
# files are installed to the proper location (provided by GNUInstallDirs)
target_include_directories(${PROJECT_NAME}
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# If we have compiler requirements for this library, list them here
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)

# Targets that we develop here
enable_testing()
add_subdirectory(test)
