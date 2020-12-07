//
//  Copyright © 2018 Arnaud Becheler. All rights reserved.
//
//
// g++ main.cpp -o main -I/usr/include/gdal -lgdal -I/home/becheler/dev -std=c++17 -I/home/becheler/dev -lboost_program_options -lsqlite3
//
//
//
#include "model_4_data.h"

int main(int argc, char* argv[])
{

  auto vm = handle_options(argc, argv);
  if (vm.count("help")) {
      return 1;
  }
  PrintVariableMap(vm);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::cout << "Initialization" << std::endl;
  SimulationContext s(vm, gen);
  std::cout << "Running ..." << std::endl;
  s.run(gen);
  return 0;
}
