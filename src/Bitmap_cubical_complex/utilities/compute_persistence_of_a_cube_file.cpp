#include <gudhi/reader_utils.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Persistent_cohomology.h>

#include "common_cube_files.h"



int main(int argc, char** argv) 
{	
	//some tests:
	//cout << "Test 1 passed? : " << test_blow_up_procedure1() << endl;
	//cout << "Test 2 passed? : " << test_blow_up_procedure2() << endl;	
	
	

  //this is a code for a single file
  std::cout
      << "This program computes persistent homology> of energy grid given in the .cube format. Please call this program with the path to the .cube file followed by three positive integers indicating how many times is the unit cell to be replicated in all directions. \n";

  if (argc != 5) {
    std::cerr << "Wrong number of parameters. Please provide the name of a file with a Perseus style bitmap at "
              << "the input. The program will now terminate.\n";
    return 1;
  }

  typedef Gudhi::cubical_complex::Bitmap_cubical_complex_periodic_boundary_conditions_base<double> Bitmap_base;
  typedef Gudhi::cubical_complex::Bitmap_cubical_complex<Bitmap_base> Bitmap_cubical_complex;
  typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
  typedef Gudhi::persistent_cohomology::Persistent_cohomology<Bitmap_cubical_complex, Field_Zp> Persistent_cohomology;
  
  unsigned blow_x = atoi( argv[4] );
  unsigned blow_y = atoi( argv[3] );
  unsigned blow_z = atoi( argv[2] );

  std::pair< std::vector<unsigned>  , std::vector< double > >  data_ = read_cube_file( argv[1] );
  
  
  std::vector<bool> directions_in_which_periodic_b_cond_are_to_be_imposed = {true, true, true};

   //std::vector<double> dataa = {0,1,2,3,4,5,6,7};
//	std::vector<unsigned> sizes = {2,2,2};
//	unsigned blow_x, blow_y, blow_z;
//	blow_x =  blow_y =  blow_z = 2;
//	std::pair< std::vector<unsigned>  , std::vector< double > >  data_ =  std::make_pair( sizes,dataa );
  
  
  std::cout << "Blowing up the unit cell. \n";
  std::pair< std::vector<unsigned>  , std::vector< double > >  data = blow_up_the_unit_cell( data_ , blow_x ,  blow_y ,  blow_z );
  cerr << "Done.\n";
  
  cout << "data.second.size() : " << data.second.size() << endl;
  
  //store the unit cell into a file:
  //store_as_cube_file( data ,  "blowed_up_cell.cube" );
  
  

  cerr << "Creating cubical complex.\n";	
  Bitmap_cubical_complex b( data.first , data.second , directions_in_which_periodic_b_cond_are_to_be_imposed );
  cerr << "Done.\n";

  // Compute the persistence diagram of the complex
  Persistent_cohomology pcoh(b);
  int p = 11;
  double min_persistence = 0;

  cerr << "Computing persistence.. \n";
  pcoh.init_coefficients(p);  // initializes the coefficient field for homology
  pcoh.compute_persistent_cohomology(min_persistence);

  std::string output_file_name(argv[1]);
  output_file_name += "_persistence";

  std::size_t last_in_path = output_file_name.find_last_of("/\\");

  if (last_in_path != std::string::npos) {
    output_file_name = output_file_name.substr(last_in_path + 1);
  }

  std::ofstream out(output_file_name.c_str());
  pcoh.output_diagram(out);
  out.close();

  std::cout << "Result in file: " << output_file_name << "\n";


//this is a code for multiple files the names of which are stored in another file:
/*
std::cout
      << "This program computes persistent homology> of energy grid given in the .cube format. Please call this program with the path to the file containg paths to files in .cube file.\n";

  if (argc != 2) {
    std::cerr << "Wrong number of parameters. Please provide the name of a file with a Perseus style bitmap at "
              << "the input. The program will now terminate.\n";
    return 1;
  }

 typedef Gudhi::cubical_complex::Bitmap_cubical_complex_periodic_boundary_conditions_base<double> Bitmap_base;
  typedef Gudhi::cubical_complex::Bitmap_cubical_complex<Bitmap_base> Bitmap_cubical_complex;
  typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
  typedef Gudhi::persistent_cohomology::Persistent_cohomology<Bitmap_cubical_complex, Field_Zp> Persistent_cohomology;
  
  
  std::vector< string > files = read_file_names_from_file( argv[1] );
  
  for ( size_t i = 0 ; i != files.size() ; ++i )
  {
	  std::pair< std::vector<unsigned>  , std::vector< double > >  data = read_cube_file( files[i].c_str() );
	  std::vector<bool> directions_in_which_periodic_b_cond_are_to_be_imposed = {true, true, true};

	  cerr << "Creating cubical complex.\n";	
	  Bitmap_cubical_complex b( data.first , data.second , directions_in_which_periodic_b_cond_are_to_be_imposed );
	  cerr << "Done.\n";

	  // Compute the persistence diagram of the complex
	  Persistent_cohomology pcoh(b);
	  int p = 11;
	  double min_persistence = 0;

	  cerr << "Computing persistence.. \n";
	  pcoh.init_coefficients(p);  // initializes the coefficient field for homology
	  pcoh.compute_persistent_cohomology(min_persistence);

	  std::string output_file_name( files[i] );
	  output_file_name += "_persistence";

	  std::size_t last_in_path = output_file_name.find_last_of("/\\");

	  if (last_in_path != std::string::npos) {
		output_file_name = output_file_name.substr(last_in_path + 1);
	  }

	  std::ofstream out(output_file_name.c_str());
	  pcoh.output_diagram(out);
	  out.close();

	  std::cout << "Result in file: " << output_file_name << "\n";
  }
*/
  return 0;
}
