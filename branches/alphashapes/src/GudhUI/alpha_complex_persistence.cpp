#include <iostream>
#include <string>


#include <QtGui/QApplication>

// to construct a Delaunay_triangulation from a OFF file
#include <gudhi/Delaunay_triangulation_off_io.h>
#include <gudhi/Alpha_complex.h>
#include <gudhi/Persistent_cohomology.h>

#include "utils/Bar_code_persistence.h"

void usage(char * const progName) {
  std::cerr << "Usage: " << progName << " filename.off " <<  // alpha_square_max_value[double] " <<
      "coeff_field_characteristic[integer > 0] min_persistence[double >= -1.0]" << std::endl;
  std::cerr << "       i.e.: " << progName << " ../../data/points/alphacomplexdoc.off 60.0 2 0.02" << std::endl;
  exit(-1); // ----- >>
}

int main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "Error: Number of arguments (" << argc << ") is not correct" << std::endl;
    usage(argv[0]);
  }
  
  QApplication qtapp(argc, argv);

  std::string off_file_name(argv[1]);
  // double alpha_square_max_value = atof(argv[2]);
  double alpha_square_max_value = 1e20;
  int coeff_field_characteristic = atoi(argv[2]);  // argv[3]
  double min_persistence = atof(argv[3]);  // argv[4]
  
  // ----------------------------------------------------------------------------
  // Init of an alpha complex from an OFF file
  // ----------------------------------------------------------------------------
  typedef CGAL::Epick_d< CGAL::Dynamic_dimension_tag > Kernel;
  Gudhi::alphacomplex::Alpha_complex<Kernel> alpha_complex_from_file(off_file_name, alpha_square_max_value);

  // ----------------------------------------------------------------------------
  // Display information about the alpha complex
  // ----------------------------------------------------------------------------
  std::cout << "Alpha complex is of dimension " << alpha_complex_from_file.dimension() <<
      " - " << alpha_complex_from_file.num_simplices() << " simplices - " <<
      alpha_complex_from_file.num_vertices() << " vertices." << std::endl;

  // Sort the simplices in the order of the filtration
  alpha_complex_from_file.initialize_filtration();

  std::cout << "Simplex_tree dim: " << alpha_complex_from_file.dimension() << std::endl;
  // Compute the persistence diagram of the complex
  Gudhi::persistent_cohomology::Persistent_cohomology< Gudhi::alphacomplex::Alpha_complex<Kernel>,
      Gudhi::persistent_cohomology::Field_Zp > pcoh(alpha_complex_from_file);
  
  std::cout << "coeff_field_characteristic " << coeff_field_characteristic <<
      " - min_persistence " << min_persistence << std::endl;

  // initializes the coefficient field for homology
  pcoh.init_coefficients(coeff_field_characteristic);

  pcoh.compute_persistent_cohomology(min_persistence);

  pcoh.output_diagram();
  
  std::vector<std::pair<double, double>> persistence_vector;
  pcoh.get_persistence(persistence_vector);

  Bar_code_persistence bc_persistence;
  
  for (auto persistence : persistence_vector) {
    bc_persistence.insert(persistence.first, persistence.second);
  }

  bc_persistence.show();

  return qtapp.exec();
}