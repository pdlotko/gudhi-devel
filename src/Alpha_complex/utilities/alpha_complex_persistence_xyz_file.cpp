
#include <boost/program_options.hpp>
#include <boost/variant.hpp>

#include <gudhi/Alpha_complex_3d.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Points_3D_off_io.h>

#include <fstream>
#include <string>
#include <vector>
#include <limits>  // for numeric_limits<>

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
* This program reads .xyz files. The first element of the returned pair is a cuboid parameter,
* the second one is the list of points.
**/ 
template< typename Point >
std::vector< Point > read_xyz_file( const char* filename )
{
	bool dbg = false;
	
	std::ifstream in;
	in.open(filename);
	if ( !in.good() )
	{
		std::cerr << "The .xyz file do not exist, program will now terminate.\n";
		throw "The .xyz file do not exist, program will now terminate.\n";
	}
	unsigned number_of_points;
	in >> number_of_points;
	
	std::vector< Point > points( number_of_points );
	
	std::string s;
	char c;
	double x,y,z;
	std::getline(in,s);
	std::getline(in,s);
	
	for ( size_t i = 0 ; i != number_of_points ; ++i )
	{
		in >> c >> x >> y >> z;
		if ( dbg )
		{
			std::cout << x << " " << y << " " << z << std::endl;
			//getchar();
		}
		points[i] = Point(x,y,z);
	}
	in.close();
	
	return points;	
}//read_xyz_file













/**
* This program reads .xyz files. The first element of the returned pair is a cuboid parameter,
* the second one is the list of points.
**/ 
std::vector< std::vector<double> > read_xyz_file_cuboid( const char* filename )
{
	std::ifstream in;
	in.open(filename);
	if ( !in.good() )
	{
		std::cerr << "The .xyz file do not exist, program will now terminate.\n";
		throw "The .xyz file do not exist, program will now terminate.\n";
	}
	unsigned number_of_points;
	in >> number_of_points;
	
	std::vector< std::vector<double> > unit_cell_coords;
	
	std::string s;	
	double x,y,z;
	std::getline(in,s);
	std::getline(in,s);
			
	//The first element in the line is a string 'Pore'. We remove it below.
	s[0] = s[1] = s[2] = s[3] = ' ';
	
	stringstream ss(s);
	ss >> x >> y >> z;
	std::vector<double> first_vect = {x,y,z};
	ss >> x >> y >> z;
	std::vector<double> second_vect = {x,y,z};
	ss >> x >> y >> z;
	std::vector<double> third_vect = {x,y,z};
	
	unit_cell_coords = {first_vect,second_vect,third_vect};
	in.close();
	
	return unit_cell_coords;	
}//read_xyz_file_cuboid




bool check_if_cuboid_is_a_cube( const std::vector< std::vector<double> >& cubo_data , double& x_min, double& y_min, double& z_min, double& x_max, double& y_max , double& z_max )
{
   //here we check if the considered cuboid is a cube:
   x_min = y_min = z_min = 0;
   x_max = y_max = z_max = 0;
  bool periodic_version = false;
  if ( (cubo_data[0][1] == cubo_data[0][2]) &&
       (cubo_data[0][2] == cubo_data[1][0]) &&
       (cubo_data[1][0] == cubo_data[1][2]) &&
       (cubo_data[1][2] == cubo_data[2][0]) &&
       (cubo_data[2][0] == cubo_data[2][1])
     )
     {
		 //in this case we know it is rectangle:
		 if ( 
		      (cubo_data[0][0] == cubo_data[1][1])&&
		      (cubo_data[1][1] == cubo_data[2][2] )
		    )
		 {
			 //and now we knot it is a cube:
			 periodic_version = true;
			 x_min = cubo_data[0][1];
			 y_min = cubo_data[1][0];
			 z_min = cubo_data[2][0];
			 x_max = cubo_data[0][0];
			 y_max = cubo_data[1][1];
			 z_max = cubo_data[2][2];
		 }
	 }
	return periodic_version;
}


int main(int argc, char **argv) 
{
  std::cout << "This program compute alpha complex persistence of a collection of points stored in .xyz file. In case the unit cell is ortogonal cube, it will impose the periodic boundary conditions. Please provide the path to .xyz file as a parameter of this program. \n\n";
	
  if ( argc != 2 )
  {
	  std::cerr << "Wrong number of parameters, the program will now terminate.\n";
	  return 1;
  }	
	
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

  
/*
  
  std::ifstream iso_cuboid_str(argv[3]);
  if (cuboid_file != std::string()) {
    if (!read_cuboid_file(cuboid_file, x_min, y_min, z_min, x_max, y_max, z_max)) {
      std::cerr << "Unable to read cuboid file " << cuboid_file << std::endl;
      exit(-1);
    }
    periodic_version = true;
  }*/
  
  

  

  

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
  double x_min, y_min, z_min, x_max, y_max , z_max;
  cout << "Will be reading points from : " << argv[1] << endl;
  std::vector< std::vector<double> > cuboid_data = read_xyz_file_cuboid( argv[1] );
  
  //cout << cuboid_data[0][0] << " "  << cuboid_data[0][1] << " " << cuboid_data[0][2] << "\n";
  //cout << cuboid_data[1][0] << " "  << cuboid_data[1][1] << " " << cuboid_data[1][2] << "\n";
  //cout << cuboid_data[2][0] << " "  << cuboid_data[2][1] << " " << cuboid_data[2][2] << "\n";
  //periodic_version = check_if_cuboid_is_a_cube( cuboid_data , x_min, y_min, z_min, x_max, y_max , z_max );
  //if ( !periodic_version )
  //{
//	  cerr << "As the bounding box is not a cuboid, the current implementation do not support periodic boundary conditions. Sorry...\n";
//  }
 // else
  //{
//	  cerr << "The bounding box is a cube, periodic bounday conditions will be imposed. \n";
//	  cout << "Here is the bounding box : " << x_min << " " << y_min << " " << z_min << " " << x_max << " " << y_max << " " << z_max << endl;
 // }
         

  switch (complexity) {
    case Gudhi::alpha_complex::complexity::FAST:
        if (periodic_version) {
          using Alpha_complex_3d =
              Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::FAST, false, true>;
            
          std::vector< Alpha_complex_3d::Point_3 > points = read_xyz_file<Alpha_complex_3d::Point_3>( argv[1] );     
          
          Alpha_complex_3d alpha_complex(points, x_min, y_min, z_min, x_max, y_max, z_max);
          alpha_complex.create_complex(simplex_tree, alpha_square_max_value);
        } else {
          using Alpha_complex_3d =
              Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::FAST, false, false>;
          std::vector< Alpha_complex_3d::Point_3 > points = read_xyz_file<Alpha_complex_3d::Point_3>( argv[1] );     
          
          Alpha_complex_3d alpha_complex(points);
          alpha_complex.create_complex(simplex_tree, alpha_square_max_value);
        }
     
      break;
    case Gudhi::alpha_complex::complexity::EXACT:
          if (periodic_version) {
          using Alpha_complex_3d =
              Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::EXACT, false, true>;
          std::vector< Alpha_complex_3d::Point_3 > points = read_xyz_file<Alpha_complex_3d::Point_3>( argv[1] );     
          
          Alpha_complex_3d alpha_complex(points, x_min, y_min, z_min, x_max, y_max, z_max);
          alpha_complex.create_complex(simplex_tree, alpha_square_max_value);
        } else {
          using Alpha_complex_3d =
              Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::EXACT, false, false>;
          std::vector< Alpha_complex_3d::Point_3 > points = read_xyz_file<Alpha_complex_3d::Point_3>( argv[1] );     
          
          Alpha_complex_3d alpha_complex(points);
          alpha_complex.create_complex(simplex_tree, alpha_square_max_value);
        }
      break;
    case Gudhi::alpha_complex::complexity::SAFE:
        if (periodic_version) {
          using Alpha_complex_3d =
              Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::SAFE, false, true>;
          std::vector< Alpha_complex_3d::Point_3 > points = read_xyz_file<Alpha_complex_3d::Point_3>( argv[1] );              
          
          Alpha_complex_3d alpha_complex(points, x_min, y_min, z_min, x_max, y_max, z_max);
          alpha_complex.create_complex(simplex_tree, alpha_square_max_value);
        } else {
          using Alpha_complex_3d =
              Gudhi::alpha_complex::Alpha_complex_3d<Gudhi::alpha_complex::complexity::SAFE, false, false>;
          std::vector< Alpha_complex_3d::Point_3 > points = read_xyz_file<Alpha_complex_3d::Point_3>( argv[1] );            
          Alpha_complex_3d alpha_complex(points);
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

	output_file_diag = std::string(argv[1]) + "_persistence";
    std::cout << "Result in file: " << output_file_diag << std::endl;
    std::ofstream out(output_file_diag);
    pcoh.output_diagram(out);
    out.close();
  
  
  

  return 0;
}
