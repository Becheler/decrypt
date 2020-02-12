cmake_minimum_required(VERSION 3.12)

find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  add_test (NAME mytest COMMAND ${Python3_EXECUTABLE} -m mytest WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()
