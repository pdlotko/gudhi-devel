#include <boost/program_options.hpp>
#include <boost/variant.hpp>

#include <gudhi/Alpha_complex_3d.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Points_3D_off_io.h>

#include <fstream>
#include <string>
#include <vector>
#include <limits>  

using namespace std;

// gudhi type definition
using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
using Filtration_value = Simplex_tree::Filtration_value;
using Persistent_cohomology =
    Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Gudhi::persistent_cohomology::Field_Zp>;



bool read_weight_file(const std::string &weight_file, std::vector<double> &weights) {
  // Read weights information from file
  std::ifstream weights_ifstr(weight_file);
  if (weights_ifstr.good()) {
    double weight = 0.0;
    // Attempt read the weight in a double format, return false if it fails
    while (weights_ifstr >> weight) {
      weights.push_back(weight);
    }
  } else {
    return false;
  }
  return true;
}

bool read_cuboid_file(const std::string &cuboid_file, double &x_min, double &y_min, double &z_min, double &x_max,
                      double &y_max, double &z_max) {
  // Read weights information from file
  std::ifstream iso_cuboid_str(cuboid_file);
  if (iso_cuboid_str.is_open()) {
    if (!(iso_cuboid_str >> x_min >> y_min >> z_min >> x_max >> y_max >> z_max)) {
      return false;
    }
  } else {
    return false;
  }
  return true;
}

template <typename AlphaComplex3d>
std::vector<typename AlphaComplex3d::Point_3> read_off(const std::string &off_file_points) {
  // Read the OFF file (input file name given as parameter) and triangulate points
  Gudhi::Points_3D_off_reader<typename AlphaComplex3d::Point_3> off_reader(off_file_points);
  // Check the read operation was correct
  if (!off_reader.is_valid()) {
    std::cerr << "Unable to read OFF file " << off_file_points << std::endl;
    exit(-1);
  }
  return off_reader.get_point_cloud();
}


/**
* The following code read out  
**/ 
template< typename Point >
std::pair< std::vector< Point > , std::vector<double> > read_xyz_atom_file( const char* filename , const std::map<string,double>& VdW )
{
	bool dbg = false;
	ifstream in;
	in.open(filename);
	
	if ( !in.good() )
	{
		cerr << "The xyz file with atom centers do not exist, program will now terminate. ";
		throw "The xyz file with atom centers do not exist, program will now terminate. ";
	}
	unsigned number_of_lines;
	in >> number_of_lines;
	std::vector< Point > points;
	std::vector<double> radius;
	points.reserve(number_of_lines);
	radius.reserve(number_of_lines);
	
	string s;
	getline(in,s);
	getline(in,s);
	string atom_type;
	double x,y,z,rad;
	
	while ( number_of_lines )
	{
		in >> atom_type >> x >> y >> z;
		
		if ( dbg )
		{
			cerr << x << " " << y << " " << z << " " << atom_type << endl;
		}
		
		--number_of_lines;
		rad = 2;
		if ( VdW.find( atom_type ) != VdW.end() )rad = VdW.find(atom_type)->second;
		
		if ( dbg )
		{
			cerr << "The corresponding VdW radius is : " << rad << endl;
			getchar();
		}
		points.push_back( Point(x,y,z) );
		radius.push_back( rad );
	}
	
	in.close();
	return std::make_pair( points , radius );
}//read_xyz_atom_file 


int main(int argc, char **argv) 
{
	
  std::string off_file_points;
  std::string weight_file;
  std::string cuboid_file;
  std::string output_file_diag;
  Filtration_value alpha_square_max_value = std::numeric_limits<Filtration_value>::infinity();
  int coeff_field_characteristic = 11;
  Filtration_value min_persistence = 0.;
  bool exact_version = false;
  bool fast_version = false;
  bool periodic_version = false;	
	
	std::cerr << "This program compute persistent homology of atom's center. We start by growing balls starting from Van der Walls balls of atoms. \n";
	std::cerr << "Please provide the name of a .xyz file as the input. \n";
	if ( argc != 2 )
	{
		cerr << "Wrong number of parameters, the program will now terminte. \n";
		throw "Wrong number of parameters, the program will now terminte. \n";
	}
	
	std::map<string,double> VdW =		
	{{"Ag",1.72},{"Ar",1.88},{"As",1.85},{"Au",1.66},{"Br",1.85},{"C", 1.70},{"Cd",1.58},{"Cl",1.75},{"Cu",1.40},{"F", 1.47},{"Ga",1.87},{"H", 1.20},{"He",1.40},{"Hg",1.55},{"I", 1.98},{"In",1.93},{"K", 2.75},{"Kr",2.02},{"Li",1.82},{"Mg",1.73},{"N", 1.55},{"Na",2.27},{"Ne",1.54},{"Ni",1.63},{"O", 1.52},{"P", 1.80},{"Pb",2.02},{"Pd",1.63},{"Pt",1.72},{"S", 1.80},{"Se",1.90},{"Si",2.10},{"Sn",2.17},{"Te",2.06},{"Tl",1.96},{"U", 1.86},{"Xe",2.16},{"Zn",1.39}};
    const char* filename = argv[1];
    
    //std::pair< std::vector< std::vector<double> > , std::vector<double> >  atoms = read_xyz_atom_file( filename , VdW );
    //std::vector<double> weights = atoms.second;
  

  double x_min = 0., y_min = 0., z_min = 0., x_max = 0., y_max = 0., z_max = 0.;
/*
  std::ifstream iso_cuboid_str(argv[3]);
  if (cuboid_file != std::string()) {
    if (!read_cuboid_file(cuboid_file, x_min, y_min, z_min, x_max, y_max, z_max)) {
      std::cerr << "Unable to read cuboid file " << cuboid_file << std::endl;
      exit(-1);
    }
    periodic_version = true;
  }
*/

  Gudhi::alpha_complex::complexity complexity = Gudhi::alpha_complex::complexity::SAFE;
  if (exact_version) {
    if (fast_version) {
      std::cerr << "You cannot set the exact and the fast version." << std::endl;
      exit(-1);
    }
    complexity = Gudhi::alpha_complex::complexity::EXACT;
  }
  if (fast_version) {
    complexity = Gudhi::alpha_complex::complexity::FAST;
  }

  Simplex_tree simplex_tree;

  switch (complexity) {
    case Gudhi::alpha_complex::complexity::FAST:
        if (periodic_version) {
          using Alpha_complex_3d =
              Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::FAST, true, true>;
              
          std::pair< std::vector< Alpha_complex_3d::Point_3 > , std::vector<double> > points_weights = 
          read_xyz_atom_file<Alpha_complex_3d::Point_3>( filename , VdW );
          
          Alpha_complex_3d alpha_complex(points_weights.first, points_weights.second, x_min, y_min, z_min, x_max, y_max, z_max);
          
          alpha_complex.create_complex(simplex_tree, alpha_square_max_value);
        } else {
          using Alpha_complex_3d =
              Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::FAST, true, false>;
          std::pair< std::vector< Alpha_complex_3d::Point_3 > , std::vector<double> > points_weights = 
          read_xyz_atom_file<Alpha_complex_3d::Point_3>( filename , VdW );
          Alpha_complex_3d alpha_complex(points_weights.first, points_weights.second);
          alpha_complex.create_complex(simplex_tree, alpha_square_max_value);
        }
      break;
    case Gudhi::alpha_complex::complexity::EXACT:
        if (periodic_version) {
          using Alpha_complex_3d =
              Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::EXACT, true, true>;
          std::pair< std::vector< Alpha_complex_3d::Point_3 > , std::vector<double> > points_weights = 
          read_xyz_atom_file<Alpha_complex_3d::Point_3>( filename , VdW );
          Alpha_complex_3d alpha_complex(points_weights.first, points_weights.second, x_min, y_min, z_min, x_max, y_max, z_max);
          alpha_complex.create_complex(simplex_tree, alpha_square_max_value);
        } else {
          using Alpha_complex_3d =
              Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::EXACT, true, false>;
          std::pair< std::vector< Alpha_complex_3d::Point_3 > , std::vector<double> > points_weights = 
          read_xyz_atom_file<Alpha_complex_3d::Point_3>( filename , VdW );
          Alpha_complex_3d alpha_complex(points_weights.first, points_weights.second);
          alpha_complex.create_complex(simplex_tree, alpha_square_max_value);
        }
      break;
    case Gudhi::alpha_complex::complexity::SAFE:
        if (periodic_version) {
          using Alpha_complex_3d =
              Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::SAFE, true, true>;
          std::pair< std::vector< Alpha_complex_3d::Point_3 > , std::vector<double> > points_weights = 
          read_xyz_atom_file<Alpha_complex_3d::Point_3>( filename , VdW );
          Alpha_complex_3d alpha_complex(points_weights.first, points_weights.second, x_min, y_min, z_min, x_max, y_max, z_max);
          alpha_complex.create_complex(simplex_tree, alpha_square_max_value);
        } else {
          using Alpha_complex_3d =
              Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::SAFE, true, false>;
          std::pair< std::vector< Alpha_complex_3d::Point_3 > , std::vector<double> > points_weights = 
          read_xyz_atom_file<Alpha_complex_3d::Point_3>( filename , VdW );
          Alpha_complex_3d alpha_complex(points_weights.first, points_weights.second);
          alpha_complex.create_complex(simplex_tree, alpha_square_max_value);
        }
      break;
    default:
      std::cerr << "Unknown complexity value " << std::endl;
      exit(-1);
      break;
  }

  // Sort the simplices in the order of the filtration
  simplex_tree.initialize_filtration();

  std::cout << "Simplex_tree dim: " << simplex_tree.dimension() << std::endl;
  // Compute the persistence diagram of the complex
  Persistent_cohomology pcoh(simplex_tree, true);
  // initializes the coefficient field for homology
  pcoh.init_coefficients(coeff_field_characteristic);

  pcoh.compute_persistent_cohomology(min_persistence);

  // Output the diagram in filediag
 
	output_file_diag = std::string(argv[1]) + "_persistence";
    std::cout << "Result in file: " << output_file_diag << std::endl;
    std::ofstream out(output_file_diag);
    pcoh.output_diagram(out);
    out.close();

  return 0;
}
