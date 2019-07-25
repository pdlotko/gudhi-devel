#include <gudhi/reader_utils.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Persistent_cohomology.h>

#include "common_cube_files.h"

std::vector< std::pair< double , double > > blow_up_diagram( const std::vector< std::pair< double , double > >& diag , unsigned bu_factor )
{
	std::vector< std::pair< double , double > > result;
	result.reserve( diag.size()*bu_factor );
	
	for ( size_t i = 0 ; i != diag.size() ; ++i )
	{
		for ( size_t j = 0 ; j != bu_factor ; ++j )
		{
			result.push_back( diag[i] );
		}
	} 
	
	return result;
}//blow_up_diagram

void store_diagram_to_file( const std::vector< std::pair< double , double > >& diag , const char* filename )
{
	ofstream out( filename );
	for ( size_t i = 0 ; i != diag.size() ; ++i )
	{
		out << diag[i].first << " " << diag[i].second << endl;
	}
	out.close();
}//store_diagram_to_file

int main(int argc, char** argv) 
{	

  //this is a code for a single file
  std::cout
      << "This program computes persistent homology of a unit cell of energy grid given in the .cube format. Subsequently the obtained persitence diagram is blwed up to obrain an approximated descriptor of the persistence homology of the blown up cell. \n"; 
      

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
  
  unsigned scaling = blow_x*blow_y*blow_z;

  std::pair< std::vector<unsigned>  , std::vector< double > >  data = read_cube_file( argv[1] );
  
  
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

  //get the persistence diagrams
  std::vector< std::pair< double , double > > diag_0 = pcoh.intervals_in_dimension( 0 );	
  std::vector< std::pair< double , double > > diag_1 = pcoh.intervals_in_dimension( 1 );	
  std::vector< std::pair< double , double > > diag_2 = pcoh.intervals_in_dimension( 2 );	
     
  //and blow them up
  std::vector< std::pair< double , double > > blowed_diag_0 = blow_up_diagram( diag_0 , scaling );	   
  std::vector< std::pair< double , double > > blowed_diag_1 = blow_up_diagram( diag_1 , scaling );	   
  std::vector< std::pair< double , double > > blowed_diag_2 = blow_up_diagram( diag_2 , scaling );
  
  //and store the results to file:
  stringstream ss0;
  ss0 << argv[1] << "_persistence_0"; 
  store_diagram_to_file( blowed_diag_0 , ss0.str().c_str() );
  
  stringstream ss1;
  ss1 << argv[1] << "_persistence_1"; 
  store_diagram_to_file( blowed_diag_1 , ss1.str().c_str() );
  
  stringstream ss2;
  ss2 << argv[1] << "_persistence_2"; 
  store_diagram_to_file( blowed_diag_2 , ss2.str().c_str() );
  	   

  return 0;
}
