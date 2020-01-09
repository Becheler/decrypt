# CMake generated Testfile for 
# Source directory: /home/becheler/dev/quetzal/test
# Build directory: /home/becheler/dev/quetzal/build/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(SpatialExpansion "model_1" "/home/becheler/dev/quetzal/test/data/europe_temp.tif")
set_tests_properties(SpatialExpansion PROPERTIES  _BACKTRACE_TRIPLES "/home/becheler/dev/quetzal/test/CMakeLists.txt;16;add_test;/home/becheler/dev/quetzal/test/CMakeLists.txt;0;")
add_test(coalescence_test "/home/becheler/dev/quetzal/test/build/testBin/coalescence_test")
set_tests_properties(coalescence_test PROPERTIES  WORKING_DIRECTORY "/home/becheler/dev/quetzal/test/build/testBin" _BACKTRACE_TRIPLES "/home/becheler/dev/quetzal/test/CMakeLists.txt;48;add_test;/home/becheler/dev/quetzal/test/CMakeLists.txt;0;")
add_test(demography_test "/home/becheler/dev/quetzal/test/build/testBin/demography_test")
set_tests_properties(demography_test PROPERTIES  WORKING_DIRECTORY "/home/becheler/dev/quetzal/test/build/testBin" _BACKTRACE_TRIPLES "/home/becheler/dev/quetzal/test/CMakeLists.txt;48;add_test;/home/becheler/dev/quetzal/test/CMakeLists.txt;0;")
add_test(expressive_test "/home/becheler/dev/quetzal/test/build/testBin/expressive_test")
set_tests_properties(expressive_test PROPERTIES  WORKING_DIRECTORY "/home/becheler/dev/quetzal/test/build/testBin" _BACKTRACE_TRIPLES "/home/becheler/dev/quetzal/test/CMakeLists.txt;48;add_test;/home/becheler/dev/quetzal/test/CMakeLists.txt;0;")
add_test(geography_test "/home/becheler/dev/quetzal/test/build/testBin/geography_test")
set_tests_properties(geography_test PROPERTIES  WORKING_DIRECTORY "/home/becheler/dev/quetzal/test/build/testBin" _BACKTRACE_TRIPLES "/home/becheler/dev/quetzal/test/CMakeLists.txt;48;add_test;/home/becheler/dev/quetzal/test/CMakeLists.txt;0;")
add_test(matrix_operation_test "/home/becheler/dev/quetzal/test/build/testBin/matrix_operation_test")
set_tests_properties(matrix_operation_test PROPERTIES  WORKING_DIRECTORY "/home/becheler/dev/quetzal/test/build/testBin" _BACKTRACE_TRIPLES "/home/becheler/dev/quetzal/test/CMakeLists.txt;48;add_test;/home/becheler/dev/quetzal/test/CMakeLists.txt;0;")
add_test(random_test "/home/becheler/dev/quetzal/test/build/testBin/random_test")
set_tests_properties(random_test PROPERTIES  WORKING_DIRECTORY "/home/becheler/dev/quetzal/test/build/testBin" _BACKTRACE_TRIPLES "/home/becheler/dev/quetzal/test/CMakeLists.txt;48;add_test;/home/becheler/dev/quetzal/test/CMakeLists.txt;0;")
