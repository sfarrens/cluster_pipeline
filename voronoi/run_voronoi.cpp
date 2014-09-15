/*
	Code:			run_voronoi.cpp
	Author:         Walter A. Santos Jr.
	Last update:	04/09/2012
	Version:		1.3
	Last changes:
		- change temporary files to stringstream
	Description: 	Code to calculate voronoi volumes and neighbors from a catalog, using voro++ routines.
	Compile:		g++ -Ivoro++/src -Wall -O3 -ansi -pedantic voronoi_walterjr.cc -o voronoi_walterjr

	Usage:	see output_help() below.
*/


// Libraries
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <cstring>

#include "wutils.h"
#include "./voro++/src/voro++.cc"

using namespace std;

// Set up the number of blocks that the container is divided
// into.
const int n_x=3,n_y=3,n_z=3; // voro++ specific parameters
const int NUM_DIMS=3; // x, y, z


void output_help(){
	cout<<"\nUsage:\n"<<flush;
	cout<<"./run_voronoi -i <xyz_file> -b <file_boundary> (-o <out_file>)\n\n"<<flush;
	cout<<"<xyz_file> is the input catalog\n"<<flush;
	cout<<"with the following columns: id x y z (...)\n\n"<<flush;
	cout<<"<file_boundary> is the file containing the overall boundaries (min/max)\n"<<flush;
	cout<<"Output file (optionally given by -o, otherwise with the tag .voro)\n"<<flush;
	cout<<"has the following columns: id x y z volume number_of_neighbors id_neighbors\n"<<flush;
}

// main
int main(int argc, char *argv[]){

	// Checking command line
		if(argc == 1){
		output_help();
		return 0;
	}
	
	// First definitions, before further checking command line
	string file_in, file_out, file_bound;
			
	for(int i=1; i<argc; i++){
		if(argv[i][0]=='-' && strlen(argv[i]) > 1){
			switch (argv[i][1]) {
				// help check
				case 'h' : 
				case 'H' :
					if(strlen(argv[i]) == 2 || (strlen(argv[i]) > 2 && string(argv[i])=="-help" ) ){
						output_help();
						return 0;
					}
					break;
				// bound file check, mandatory option
				case 'b' :
					if(strlen(argv[i]) == 2 || (strlen(argv[i]) > 2 && string(argv[i])=="-bound" ) ){
						if((i+1) < argc && argv[i+1][0]!='-')
							file_bound = argv[i+1];
						if((i+2) < argc && argv[i+2][0]!='-'){
							cout << "ERROR: Only one boundary file is accpeted.\n\n";
							output_help();
							exit(-1);
						}
					}
					break;
				// input xyz file check, mandatory option
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
					break;
				// input out file, optional
				case 'o' :
					if((i+1) < argc && strlen(argv[i]) == 2 || (strlen(argv[i]) > 2 && string(argv[i])=="-out" ) ){
						if((i+1) < argc && argv[i+1][0]!='-')
							file_out = argv[i+1];
						if((i+2) < argc && argv[i+2][0]!='-'){
							cout << "ERROR: Only one output file is accepted.\n\n";
							output_help();
							exit(-1);
						}
					}
					break;
			}
			
		}
	}
	
	// checking for boundary file, mandatory
	if(file_bound.length() == 0){
		cout << "ERROR: Input boundary file is mandatory.\n\n";
		output_help();
		exit(-1);
	}
	
	// checking for xyz file, mandatory
	if(file_in.length() == 0){
		cout << "ERROR: Input xyz file is mandatory.\n\n";
		output_help();
		exit(-1);
	}
	
	// reading and saving boundary data
	vector<double> min, max;
	min.resize(NUM_DIMS);
	max.resize(NUM_DIMS);
	int _cols, _rows;
	vector< vector<string> > bound_mx;
	
	wutils::read_ascii(file_bound, bound_mx, _rows, _cols);
	for(int k=0; k<NUM_DIMS; k++){
		min[k] = atof(bound_mx[0][k].c_str());
		max[k] = atof(bound_mx[1][k].c_str());
	}
	
	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block.
	container con(min[0],max[0],min[1],max[1],min[2],max[2],n_x,n_y,n_z,
			false,false,false,8);
			
	// reading and saving xyz data, it needs a temp file with just id x y z
	int ncols, nrows;
	vector< vector<string> > data_mx;
	
	//wutils::read_ascii(file_in, data_mx, nrows, ncols);
	wutils::readFile(file_in, data_mx, nrows, ncols);
	//string file_in_temp = wutils::tag_outfile(file_in, "temp");
	stringstream in_temp (stringstream::in | stringstream::out);
	
	int id = 0; // Have to create a new int id. :/
	// saving id x y z into a new temp file
	for (int i=0; i<nrows; i++){
		in_temp << id << " ";
		for(int j=1; j<NUM_DIMS+1; j++)
			in_temp << data_mx[j][i] << " ";
		in_temp << "\n";
		id++;
	}
	//in_temp.close();
	in_temp.seekg (ios::beg);
	
	
	// string test;
	// for(int i=0;i<5;i++){
	// in_temp >> test;
	// cout << test << endl;
	// }
	
	// Import the monodisperse test packing and output the Voronoi
	// tessellation in gnuplot and POV-Ray formats.
	con.import(in_temp);
			
	// Now it also needs a temp out file
	//string file_out_temp = wutils::tag_outfile(file_in, "out.temp");
	stringstream out_temp (stringstream::in | stringstream::out);
	
	// Do a custom output routine that outputs the particle IDs and
	// positions, plus the volume, number of neighbors and neighbors IDS
	con.print_all_custom("%i %q %v %s %n", out_temp);
	
	out_temp.seekg (ios::beg);
	// string test;
	// for(int i=0;i<5;i++){
	// out_temp >> test;
	// cout << test << endl;
	 // }
	
	
	// now join the output with the input info that was left aside,
	// removing the temp files
	
	// first read the voro++ output
	//ifstream out_temp (file_out_temp.c_str());
	
	ofstream out;
	if(file_out.length() > 0)
		out.open(file_out.c_str());
	
	int num_neighbors;
	double aux, volume;
		
	for(int i=0; i < nrows; i++){
		//read voro intrinsic id, unimportant aux, volume and num neighbors
		out_temp >> id;
		for(int j=1; j<4; j++)
			out_temp >> aux;
		out_temp >> volume;
		out_temp >> num_neighbors;
		
		// if out to a file...
		if(file_out.length() > 0){
			out << data_mx[0][id] << " ";
			for(int j=1; j<=5; j++)
				out << data_mx[j][id] << " ";
			out << 1/volume << " "; //density
			out << num_neighbors << " ";
			// saving neighbor ids
			for(int j=6; j <= 5+num_neighbors; j++){
				out_temp >> aux;
				out << aux << " ";
			}
			out << "\n";
		}
		else { // same as above, but to cout
			cout << data_mx[0][id] << " ";
			for(int j=1; j<=5; j++)
				cout << data_mx[j][id] << " ";
			cout << 1/volume << " "; //density
			cout << num_neighbors << " ";
			// saving neighbor ids
			for(int j=6; j <= 5+num_neighbors; j++){
				out_temp >> aux;
				cout << aux << " ";
			}
			cout << "\n";
		}
	}
	//out_temp.close();
	if(file_out.length() > 0){
		cout << "Output saved in " << file_out << "." << endl;
		out.close();
	}
	
	//remove(file_in_temp.c_str());
	//remove(file_out_temp.c_str());
	
	return 0;
}
