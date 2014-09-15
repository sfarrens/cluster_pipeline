/*
	Code:			run_splitradec.cpp
	Author:         Walter A. Santos Jr. (walter.augusto@gmail.com)
					(Based on code from Fotini Economou)
	Last update:	09/09/2013
	Version: 		1.3
	Last changes: 	- Bug fixes with 'neighbors' healpix function.
	Description: 	Code to split a radec space of galaxies into smaller intersecting regions.
	Usage:			see output_help() below
	Compile (in linux): g++ -o run_splitradec run_splitradec.cpp wutils.cpp -I/home/cosmos/library_src/Healpix_2.20a/src/cxx/linux_icc/include -L/home/cosmos/library_src/Healpix_2.20a/src/cxx/linux_icc/lib -L/state/partition1/apps/intel/Compiler/11.1/046/lib/intel64/ -lhealpix_cxx -lcxxsupport -lpsht -lc_utils -lfftpack -fopenmp -lirc -limf					
*/

#include <iostream>
#include <fstream>
#include <vector>
//#include <arr.h>

#include "/home/cosmos/library/include/Healpix_2.20a/xcomplex.h"
#include "/home/cosmos/library/include/Healpix_2.20a/cxxutils.h"
#include "/home/cosmos/library/include/Healpix_2.20a/paramfile.h"
// #include "../healpix/cxxsupport/simparams.h"
#include "/home/cosmos/library/include/Healpix_2.20a/healpix_data_io.h"
#include "/home/cosmos/library/include/Healpix_2.20a/alm.h"
#include "/home/cosmos/library/include/Healpix_2.20a/healpix_map.h"
#include "/home/cosmos/library/include/Healpix_2.20a/healpix_map_fitsio.h"
// #include "/home/cosmos/library/include/Healpix_2.20a/alm_map_tools.h"
#include "/home/cosmos/library/include/Healpix_2.20a/alm_powspec_tools.h"
#include "/home/cosmos/library/include/Healpix_2.20a/fitshandle.h"
#include "/home/cosmos/library/include/Healpix_2.20a/healpix_base.h"

#include "wutils.h"

using namespace std;

const Healpix_Ordering_Scheme SCHEME = NEST;
//const double DEFAULT_INTERS = 0.10; //percentage
const int DEFAULT_NUM_PIXELS = 3072; // default number of pixels that the *whole sky* should be divided into
const int DEFAULT_NUM_PIXELS_INTERS = 49152;

const int array_pixels[] = {12,48,192,768,3072,12288,49152,196608,786432,3145728,12582912};
const vector<int> pixels_values (array_pixels, array_pixels + sizeof(array_pixels) / sizeof(int) );

// TODO: UPDATE OUTPUT_HELP
void output_help(){
	cout<<"\nUsage:\n"<<flush;
	cout<<"./run_splitradec -i <radec_file> (-p <num_pixels>) (-inters <perc>) (-o <file_split_n>) -n <n>\n\n"<<flush;
	cout<<"where <num_pixels> is the number of pixels that the *whole sky* should be divided into,\n\n"<<flush;
	cout<<"<radec_file> is the input catalog\n"<<flush;
	cout<<"with the following columns: id ra dec (...)\n\n"<<flush;
	cout<<"n is the output file int number to split to (0-<num_pixels>).\n"<<flush;
	cout<<"<file_split_n> is the 'split' output file associated with n. It probably must be the same name\n"<<flush;
	cout<<"if you run for others <radec_file> with the same n. If it's not given, it outputs to cout.\n\n"<<flush;
	cout<<"<perc> is the percentage of intersection between sub-spaces (optional)\n\n"<<flush;
}

/*
12
48
192
768
3072
12288
49152
196608
786432
3145728
12582912
*/

int  main(int argc, char *argv[]){

	// Checking command line
	if(argc == 1){
		output_help();
		return 0;
	}
	
	int num_split_pixs=-1, num_inters_pixs=-1;
	double inters=-1;
	string file_in;
	string file_out;
	//string file_bounds, file_out;
	int outn = -1;
			
	for(int i=1; i<argc; i++){
		if(argv[i][0]=='-' && strlen(argv[i]) > 1){
			switch (argv[i][1]) {
				case 'h' : 
				case 'H' :
					if(strlen(argv[i]) == 2 || (strlen(argv[i]) > 2 && string(argv[i])=="-help" ) ){
						output_help();
						return 0;
					}
					break;
				case 'i' :
					if(strlen(argv[i]) == 2 || (strlen(argv[i]) > 2 && string(argv[i])=="-inp" ) ){
						
						if((i+1) < argc && argv[i+i][0]!='-'){
							file_in = argv[i+1];
							
							}
						if((i+2) < argc && argv[i+2][0]!='-'){
							cout << "ERROR: Only one input file is accepted.\n\n";
							output_help();
							exit(-1);
						}
					}
					else if(strlen(argv[i]) > 2 && string(argv[i])=="-inters"){
							if(argv[i+1][0]!='-' && atof(argv[i+1]) >= 0 && atof(argv[i+1]) < 1){
								inters = atof(argv[i+1]);
								
							}
					}
					break;
				case 'o' :
					if((i+1) < argc && (strlen(argv[i]) == 2 || (strlen(argv[i]) > 2 && string(argv[i])=="-out"))){
						// checking for output file
						file_out = argv[i+1];
						if((i+2) < argc && argv[i+2][0]!='-'){
							cout << "ERROR: Only one output file is accepted.\n\n";
							output_help();
							exit(-1);
						}
					}
					break;
				case 'n' :
					if(strlen(argv[i]) == 2 ){
						if((i+1) < argc && argv[i+1][0]!='-'){
							outn = atoi(argv[i+1]);
						}
					}
					break;
				case 'p' :
					if(strlen(argv[i]) == 2){
						if((i+1) < argc && argv[i+1][0]!='-' && atoi(argv[i+1]) > 0){
							num_split_pixs = atoi(argv[i+1]);
						}
					}
					break;					
			}
		}
	}

	if(file_in.length() == 0){
		cout << "ERROR: Input file is mandatory.\n\n";
		output_help();
		exit(-1);
	}
	
	if(num_split_pixs <= 0)
		num_split_pixs = DEFAULT_NUM_PIXELS;
	else {
		int i;
		for(i=0; i<pixels_values.size(); i++)
			if(num_split_pixs <= pixels_values[i]){
				num_split_pixs = pixels_values[i];
				break;
			}
		if(i==pixels_values.size())
			num_split_pixs = pixels_values.back();
	}
	
	if(inters < 0 || inters > 1)
		num_inters_pixs = DEFAULT_NUM_PIXELS_INTERS;	
	else{
		num_inters_pixs = (int) ceil(num_split_pixs/(inters*inters));
		int i;
		for(i=0; i<pixels_values.size(); i++)
			if(num_inters_pixs <= pixels_values[i]){
				num_inters_pixs = pixels_values[i];
				break;
			}
		if(i==pixels_values.size())
			num_inters_pixs = pixels_values.back();
	}
	
	/*
	int _cols, _rows;
	vector< vector<string> > in_mx;
	wutils::read_ascii(file_bounds, in_mx, _rows, _cols);
	
	string string_split = "Number splitting pixels:";
	string string_inters = "Number intersecting pixels:";
	for(int i=0; i<_rows; i++){
		if (in_mx[1][i].compare("splitting") == 0)
			num_split_pixs = atoi(in_mx[3][i].c_str());
		if (in_mx[1][i].compare("intersecting") == 0)
			num_inters_pixs = atoi(in_mx[3][i].c_str());
	}
	*/
	
	if(outn < 0 || outn > num_split_pixs-1){
		cout << "ERROR: Output file number was not found or its input value was wrong.\n\n";
		output_help();
		exit(-1);
	}
	
	int nside_split = (int) sqrt(num_split_pixs/12);
	int nside_inters = (int) sqrt(num_inters_pixs/12);
		
	class Healpix_Base heal_split_pix;
	class Healpix_Base heal_inters_pix;
	heal_split_pix.SetNside(nside_split, SCHEME);
	heal_inters_pix.SetNside(nside_inters, SCHEME);
	
	//cout << outn << " " << heal_split_pix.pix2ang(outn) << endl;
	
	// getting the small inters pixels *inside* the outn split pix
	// and also getting the inters region using those *inside* pix
	vector<int> vec_inters_pix_inside;
	vector<int> vec_inters_region_pix;
	fix_arr< int, 8 > neighbors_inters_pix;
	for(int i=0; i<num_inters_pixs; i++){
		class pointing inters_pix_ang;
		inters_pix_ang = heal_inters_pix.pix2ang(i);
		int inters_pix_as_split = heal_split_pix.ang2pix(inters_pix_ang);
		if (inters_pix_as_split == outn){
			vec_inters_pix_inside.push_back(i);
			heal_inters_pix.neighbors(i, neighbors_inters_pix);
			// now see if the neighbors are inside or not. If not, they are in the inters region
			for(int j=0; j<8; j++){
				if(neighbors_inters_pix[j]>-1){
					class pointing inters_neigh_pix_ang = heal_inters_pix.pix2ang(neighbors_inters_pix[j]);
					int inters_pix_neigh_as_split = heal_split_pix.ang2pix(inters_neigh_pix_ang);
					if(inters_pix_neigh_as_split != outn)
						vec_inters_region_pix.push_back(neighbors_inters_pix[j]);
				}
			}
		}
	}
	
	int ncols, nrows;
	vector< vector<string> > data_mx;
	wutils::read_ascii(file_in, data_mx, nrows, ncols);
	
	double ra, dec;
	class pointing gal_coord;
		
	//open out file
	ofstream out;
	if(file_out.length() > 0)
		out.open(file_out.c_str(), ios::out | ios::app);
  
	for(int i=0; i<nrows; i++){
	
		//vector<int> listpix;

		ra = atof(data_mx[1][i].c_str());		// TODO: check later columns for RA DEC in input file
		dec = atof(data_mx[2][i].c_str());
		
		ra*=M_PI/180.0;
		dec*=M_PI/180.0;
    
		gal_coord.theta = M_PI/2. - dec ; //theta=90-dec
		gal_coord.phi=ra;  
		
		// now see if the point is in either the inside or neighboring inters pixels
		int gal_inters_pix = heal_inters_pix.ang2pix(gal_coord);
		bool isinside = false;
		for (int j=0; j<vec_inters_pix_inside.size(); j++){
			if(gal_inters_pix == vec_inters_pix_inside[j]){
				isinside = true;
				break;
			}
		}
		if (!isinside){
			for (int j=0; j<vec_inters_region_pix.size(); j++){
				if(gal_inters_pix == vec_inters_region_pix[j]){
					isinside = true;
					break;
				}
			}
		}
		
		if(isinside){
			if(file_out.length() > 0){
				for(int k=0; k<ncols; k++)
					out << data_mx[k][i] << " ";
				out << "\n";
				//out << file_in << " " << i << "\n";
			}
			else{
				for(int k=0; k<ncols; k++)
					cout << data_mx[k][i] << " ";
				cout << "\n";
				//cout << file_in << " " << i << "\n";
			}
		}
		
	}
	
	if(file_out.length() > 0)
		out.close();

	return 0;
}