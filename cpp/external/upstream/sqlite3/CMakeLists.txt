cmake_minimum_required(VERSION 3.12)

find_package(SQLite3 CONFIG QUIET)

if(SQLite3_FOUND)
  get_property(_loc TARGET SQLite3::sqlite3 PROPERTY LOCATION)
  message(STATUS "Found SQLite3: ${_loc} (found version ${SQLite3_VERSION})")
  add_library(sqlite3_external INTERFACE)  # dummy
else()
  message(STATUS "Suitable SQLite3 could not be located. Downloading and building")

  # Add CMake project file
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/temp)
  file(WRITE ${CMAKE_BINARY_DIR}/temp/CMakeLists.txt
      "cmake_minimum_required(VERSION 3.12)\n"
      "set(PROJECT_NAME sqlite3_external)\n"
      "include_directories(${CMAKE_BINARY_DIR}/source)\n"
      "add_library(sqlite3 ${CMAKE_BINARY_DIR}/source/sqlite3.c)\n"
      "install(TARGETS sqlite3 DESTINATION lib)\n"
      "install(FILES sqlite3.h DESTINATION include)\n")

  # Download, configure, build and install.
  include(ExternalProject)
  ExternalProject_Add(sqlite3_external
      URL               https://www.sqlite.org/2020/sqlite-autoconf-3310100.tar.gz
      URL_HASH          SHA1=0c30f5b22152a8166aa3bebb0f4bc1f3e9cc508b
      DOWNLOAD_NO_PROGRESS  1
      #--Update/Patch step----------
      # UPDATE_COMMAND    ""
      UPDATE_COMMAND    ${CMAKE_COMMAND} -E copy
                        ${CMAKE_BINARY_DIR}/temp/CMakeLists.txt
                        <SOURCE_DIR>/CMakeLists.txt
      #--Configure step-------------
      #--Build step-----------------
      BUILD_COMMAND     ${CMAKE_COMMAND} --build .
      LOG_BUILD         1
      BUILD_IN_SOURCE   1
    )

endif()
