// standard stuff
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstddef>



using namespace std;


/**
* This procedure takes a counter of a length 3 as well as the final counter, and increment the first counter if only possible.
* If the first counter gets incremented, the procedure do it and return true. In the other case it does nothing and return false.
**/
bool increment_counter( std::vector<int>& counter_to_increment , const std::vector<int>& end_counter )
{
	if ( counter_to_increment[0] < end_counter[0]-1 )
	{
		counter_to_increment[0]++;
	}
	else
	{
		counter_to_increment[0] = 0;
		if ( counter_to_increment[1] < end_counter[1]-1 )
		{
			counter_to_increment[1]++;
		}
		else
		{
			counter_to_increment[1] = 0;
			if ( counter_to_increment[2] < end_counter[2]-1 )
			{
				counter_to_increment[2]++;
			}
			else
			{
				return false;
			}
		}
	}
	return true;
}//increment_counter


/**
 * Given the current and the end counter this procedure return the position in the vector with the corresponding value.
**/
unsigned find_position_corresponding_to_counter( const vector<int>& counter , const std::vector<int>& end_counter  )
{
	return counter[0] + end_counter[0]*counter[1]  + end_counter[0]*end_counter[1]*counter[2];
}//find_position_corresponding_to_counter


/**
 * Given the current and the end counter this procedure return the position in the vector with the corresponding value.
**/
unsigned find_position_corresponding_to_counter( const vector<unsigned>& counter , const std::vector<unsigned>& end_counter  )
{
	return counter[0] + end_counter[0]*counter[1]  + end_counter[0]*end_counter[1]*counter[2];
}//find_position_corresponding_to_counter

std::vector< int > find_counter_corresponding_to_position( unsigned position , std::vector< int > params )
{
	std::vector< int >  result(3);
	
	cerr << "position : "<<position << "(params[0]*params[1]) : " << (params[0]*params[1]) << endl;
	
	result[2] = position/(params[0]*params[1]);
	position = position%(params[0]*params[1]);
	
	
	
	result[1] = position/params[0];
	position =  position%params[0];
	
	result[0] = position;
	
	
	return result;
}//find_counter_corresponding_to_position



void output_counter( const std::vector<int>& counter )
{
	cout << counter[0] << " " << counter[1] << " " << counter[2] << endl;
}//output_counter



std::pair< std::vector<unsigned>  , std::vector< double > > blow_up_the_unit_cell( const std::pair< std::vector<unsigned>  , std::vector< double > > uc , unsigned blow_x ,  unsigned blow_y ,  unsigned blow_z )
{
	bool dbg = false;
	std::vector<int> current = {0,0,0};
	//std::vector<int> end = {blow_x*(uc.first[0]+1)-1 , blow_y*(uc.first[1]+1)-1 , blow_z*(uc.first[2]+1)-1 };
	std::vector<int> end = {(int)(blow_x*uc.first[0]) , (int)(blow_y*uc.first[1]) , (int)(blow_z*uc.first[2]) };

	//if ( dbg )
	{
		cerr << "We will blow up the cells from : " << uc.first[0] << " " << uc.first[1] << " " << uc.first[2] << "  to the size : " << end[0] << " , " << end[1] << " , " << end[2] << endl;
		//getchar();
	}

	std::vector< double > blowed_up_data;
	blowed_up_data.reserve( end[0]*end[1]*end[2] );
	std::vector< unsigned > modulo_counter = {0,0,0};

	while ( true )
	{
		modulo_counter[0] = current[0]%uc.first[0];
		modulo_counter[1] = current[1]%uc.first[1];
		modulo_counter[2] = current[2]%uc.first[2];

		if ( dbg )
		{
			output_counter(current);
			cout << "Cunter's position : " << find_position_corresponding_to_counter(current, end) << endl;
			cout << "Modulo counter : " << modulo_counter[0] << " " << modulo_counter[1] << " " << modulo_counter[2] << endl;
			cout << "Position corresponding to modulo counter : " << find_position_corresponding_to_counter( modulo_counter , uc.first ) << endl;
			getchar();
		}

		blowed_up_data.push_back( uc.second[ find_position_corresponding_to_counter( modulo_counter , uc.first ) ] );

		if ( !increment_counter( current , end ) )break;
	}


	std::vector<unsigned> endd = {blow_x*uc.first[0] , blow_y*uc.first[1] , blow_z*uc.first[2] };

	return std::make_pair( endd , blowed_up_data );
}//blow_up_the_unit_cell



/**
 * This function read .cube file and return a pair consisitng of the sizes of the
 * grid and vector with all filtration values.
**/
std::pair< std::vector<unsigned>  , std::vector< double > > read_cube_file( const char* filename )
{
	ifstream in;
	in.open( filename );

	if ( !in.good() )
	{
		cerr << "The file you have provided do not exist. The program will now terminate. \n";
		throw "The file you have provided do not exist. The program will now terminate. \n";
	}
	
	//read the first two lines
	string s;
	getline(in,s);
	getline(in,s);

	double buffer;
	unsigned x_size, y_size, z_size;

	//first we read the number of atoms:
	unsigned number_of_atoms;
	in >> number_of_atoms >> buffer >> buffer >> buffer >> x_size >> buffer >> buffer >> buffer >> y_size >> buffer >> buffer >> buffer >> z_size;
	cerr << "We have : " << number_of_atoms << " atoms. The size of the grid is : " << x_size << " , " << y_size << " , " << z_size << endl;

	std::vector< unsigned > params_of_grid = { z_size , y_size , x_size };
	unsigned size_of_grid = x_size * y_size * z_size;

	//now we need to skip number_of_atoms lines
	string line;
	for ( size_t i = 0 ; i <= number_of_atoms ; ++i )
	{
		std::getline( in,line );
	}

	std::vector< double > filtration;
	filtration.reserve(size_of_grid);
	double value;
	while ( size_of_grid )
	{
		in >> value;
		filtration.push_back( value );
		--size_of_grid;
	}

	return std::make_pair( params_of_grid , filtration );
}//read_cube_file


/**
 * Suppose that we have a large number of files in a folder, and we want to run the code for each of them. We can get the list of those files with ls command and put its result to a file.
 * This procedure read such a file and return vector of files, so that the code can be dun for each of them in a loop.
**/
std::vector< string > read_file_names_from_file( const char* filename )
{
	ifstream in;
	in.open( filename );
	std::vector< string > result;

    typedef string::size_type string_size;


    // invariant: we have processed characters [original value of i, i)
    while ( in.good() )
    {
		string_size i = 0;
        string s;
		getline( in,s );
		while (i != s.size())
		{
		   // ignore leading blanks
		   // invariant: characters in range [original i, current i) are all spaces
		   while (i != s.size() && isspace(s[i]))++i;

		   // find end of next word
		   string_size j = i;
		   // invariant: none of the characters in range [original j, current j)is a space
		   while (j != s.size() && !isspace(s[j]))j++;

		   // if we found some nonwhitespace characters
		   if (i != j) {
			  // copy from s starting at i and taking j - i chars
			  result.push_back(s.substr(i, j - i));
			  i = j;
		   }
	   }
   }
	in.close();
	return result;
}

bool test_blow_up_procedure1()
{
	std::vector<double> data = {0,1,2,3,4,5,6,7};
	std::vector<unsigned> sizes = {2,2,2};

	std::vector< double > template_solution = {0,1,0,1,2,3,2,3, 0,1,0,1,2,3,2,3,  4,5,4,5,6,7,6,7, 4,5,4,5,6,7,6,7, 0,1,0,1,2,3,2,3, 0,1,0,1,2,3,2,3,  4,5,4,5,6,7,6,7, 4,5,4,5,6,7,6,7};

	std::pair< std::vector<unsigned>  , std::vector< double > > bj = blow_up_the_unit_cell( std::make_pair( sizes,data ) , 2,2,2 );


	for ( size_t i = 0 ; i != bj.second.size() ; ++i )
	{
		//cout << bj.second[i] << " ";
		//if ( (i+1)%(bj.first.size()+1)==0 )cout << endl;
		if ( bj.second[i] != template_solution[i] ) return false;
	}
	return true;
}//test_blow_up_procedure1



bool test_blow_up_procedure2()
{
	std::vector<double> data = {0,1,2,3,4,5,6,7,8, 9,10,11,12,13,14,15,16,17, 18,19,20,21,22,23,24,25,26};
	std::vector<unsigned> sizes = {3,3,3};

	std::vector< double > template_solution = {0,1,2,0,1,2,3,4,5,3,4,5,6,7,8,6,7,8,0,1,2,0,1,2,3,4,5,3,4,5,6,7,8,6,7,8,9,10,11,9,10,11,12,13,14,12,13,14,15,16,17,15,16,17,9,10,11,9,10,11,12,13,14,12,13,14,15,16,17,15,16,17,18,19,20,18,19,20,21,22,23,21,22,23,24,25,26,24,25,26,18,19,20,18,19,20,21,22,23,21,22,23,24,25,26,24,25,26,0,1,2,0,1,2,3,4,5,3,4,5,6,7,8,6,7,8,0,1,2,0,1,2,3,4,5,3,4,5,6,7,8,6,7,8,9,10,11,9,10,11,12,13,14,12,13,14,15,16,17,15,16,17,9,10,11,9,10,11,12,13,14,12,13,14,15,16,17,15,16,17,18,19,20,18,19,20,21,22,23,21,22,23,24,25,26,24,25,26,18,19,20,18,19,20,21,22,23,21,22,23,24,25,26,24,25,26};

	std::pair< std::vector<unsigned>  , std::vector< double > > bj = blow_up_the_unit_cell( std::make_pair( sizes,data ) , 2,2,2 );


	for ( size_t i = 0 ; i != bj.second.size() ; ++i )
	{
		//cout << bj.second[i] << " ";
		//if ( (i+1)%6==0 )cout << endl;
		if ( bj.second[i] != template_solution[i] ) return false;
	}
	return true;
}//test_blow_up_procedure2

/**
* This procedure simply open a .cube file and read the first number, which is the number of atoms
**/
unsigned read_number_of_atoms_from_cube_file( const char* filename )
{
	ifstream in;
	in.open( filename );

	if ( !in.good() )
	{
		cerr << "The file you have provided do not exist. The program will now terminate. \n";
		throw "The file you have provided do not exist. The program will now terminate. \n";
	}

	//first we read the number of atoms:
	unsigned number_of_atoms;
	in >> number_of_atoms;
	in.close();

	return number_of_atoms;
}//read_number_of_atoms_from_cube_file

/**
* This procedure stores the given data as a file in a .cube format.
**/
void store_as_cube_file( const std::pair< std::vector<unsigned>  , std::vector< double > >& data ,  const char* filename )
{
	ofstream out( filename );
	out << "Blowed up cell" << endl;
	out << "Blowed up cell" << endl;
	out << "1 0 0 0" << endl;
	out << data.first[2] << " 1 0 0" << endl;
	out << data.first[1] << " 0 1 0" << endl;
	out << data.first[0] << " 0 0 1" << endl;
	out << "1 0 0 0 0" << endl;

	size_t num = 0;
	for ( size_t i = 0 ; i != data.second.size() ; ++i )
	{
		++num;
		out << data.second[i];
		if ( i != data.second.size()-1 )
		{
			out << " ";
		}
		//if ( num%6 ==0 )out << endl;
		if ( num == data.first[0] )
		{
			out << endl;
			num = 0;
		}
	}

	out.close();
}//store_as_cube_file

