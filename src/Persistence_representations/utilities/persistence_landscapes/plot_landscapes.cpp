/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Persistence_landscape.h>

#include <iostream>
#include <sstream>

using Persistence_landscape = Gudhi::Persistence_representations::Persistence_landscape;

int main(int argc, char** argv) {
  std::cout << "This program creates a gnuplot script from a persistence landscape stored in a file (the file needs "
            << "to be created beforehand). Please call the code with the name of a single landscape file.\n";
  if (argc != 2) {
    std::cout << "Wrong parameter list, the program will now terminate \n";
    return 1;
  }

  Persistence_landscape l;
  l.load_landscape_from_file(argv[1]);
  l.plot(argv[1]);

  return 0;
}
