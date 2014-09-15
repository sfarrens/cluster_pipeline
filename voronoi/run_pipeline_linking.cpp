/*
	Code:			run_pipeline_linking.cpp
	Author:         Walter A. Santos Jr. (walter.augusto@gmail.com)
	Last update:	08/08/2013
	Version: 		1.3
	Changes: 		- Bug fixes
					- Fixing output_help()
					- Includes several new comments
					- new input parameter for the pre and post path to the files in the given filename
	TODO: 			Testing in the pipeline; try to not assume id/ra/dec column positions.
	Description: 	Code to link voronoi3D results to the next steps (sam's part) in the pipeline,
					by tracing back radec.
	Usage:			see output_help() function below.
*/

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

// #include "wutils.h"

using namespace std;

void output_help(){
	cout<<"\nUsage:\n"<<flush;
	cout<<"./run_pipeline_linking -inp <file> -f <filenames> (-pre <pre_path>) (-post <post_path>) (-o <file_out>)\n\n"<<flush;
	cout<<"where the input <file> is the output of the voro out results\n"<<flush;
	cout<<"with the following columns: id x y z orig_file orig_line density num_neighbors ids_neighbors(...)\n\n"<<flush;
	cout<<"<filenames> is an ascii file that has the original input filenames (read from params ini file).\n"<<flush;
	cout<<"<path> is the base path that should be applied to the filename(s). Optional if not needed.\n\n"<<flush;
	cout<<"A new (optional, file_out) file is generated with: id ra dec redshift density.\n\n"<<flush;
}

//counting the rows / total numer of points
int count_num_points (std::string filename){
	int size=0;
	std::string line;

	std::ifstream file (filename.c_str());
	
	if (!(file.is_open())){
		std::cout << "ERROR: cannot open file "<< filename << std::endl;
		exit(-1);
	}

	while( getline(file, line) )
		size++;
   
	file.close();
	return size;
}

// function to trace back the data (id/ra/dec/redshift), given originals filename and row number
void lookbackCat(std::string file_voro, std::string filenames, std::string pre_path, std::string post_path, std::string file_out){
	int num_files = count_num_points(filenames);
	int num_total = count_num_points(file_voro);
	
	std::ifstream v (file_voro.c_str());
	std::ifstream fls (filenames.c_str());
	
	if (!(v.is_open()) || !(fls.is_open())){
		std::cout << "ERROR: cannot open files "<< file_voro << " and/or " << filenames << "to look back!" << std::endl;
		exit(-1);
	}
	
	ofstream out;
	if(file_out.length() > 0)
		out.open(file_out.c_str());
	
	std::ifstream fk;
	std::string filenamek, filename_in, id, line;
	int nrow, nmatches=0;
	double aux, density, ra, dec, redshift;
	
	// looping through original filenames
	for (int k=0; k<num_files; k++){
		fls >> filenamek;
		filenamek = pre_path + filenamek + post_path;
		fk.open(filenamek.c_str());
		
		if(fk.is_open()){
			for(int i=0; i<num_total; i++){
				v >> id;
				for(int j=0; j<3; j++)
					v >> aux;
				v >> filename_in;
				// if we have a match in the filename, read row number,
				// read nrows to get to the specific line, then
				// save radec and redshift
				if (filename_in.compare(filenamek) == 0){ //matched filenames
					nmatches++;
					v >> nrow;
					v >> density;
					for(int ki=0; ki<nrow; ki++)
						getline (fk, line);
					fk >> aux;
					fk >> ra;
					fk >> dec;
					fk >> redshift;
					if(file_out.length() > 0)
						out << id << " " << ra << " " << dec << " " << redshift << " " << density << "\n";
					cout << id << " " << ra << " " << dec << " " << redshift << " " << density << "\n";
				}
				getline (v, line);
				fk.seekg (0, std::ios::beg);
			}
		fk.close();
		
		}
		else{
			std:cout << "WARNING: cannot open file " << filenamek << "." << std::endl;
		}
		// just time saving condition... if all the points from voro input are read,
		// then stop the original filenames loop
		if(nmatches >= num_total)
			break;
		v.seekg (0, std::ios::beg);
	}
	
	v.close();
	fls.close();
	
	if(file_out.length() > 0){
		cout << "Output saved in " << file_out << "." << endl;
		out.close();
	}

}


int main(int argc, char *argv[]){
	string file_in, file_filenames, file_out, pre_path, post_path;
		
	// Checking command line
	if(argc == 1){
		output_help();
		return 0;
	}
	for(int i=1; i<argc; i++){
		if(argv[i][0]=='-' && strlen(argv[i]) > 1){
			switch (argv[i][1]) {
				// checking for help option
				case 'h' : 
				case 'H' :
					if(strlen(argv[i]) == 2 || (strlen(argv[i]) > 2 && string(argv[i])=="-help" ) ){
						output_help();
						return 0;
					}
					break;
				// checking for input (voro) file	
				case 'i' :
					if(strlen(argv[i]) == 2 || (strlen(argv[i]) > 2 && string(argv[i])=="-inp" ) ){
						if((i+1) < argc && argv[i+1][0]!='-'){
							file_in = argv[i+1];
						}
						if((i+2) < argc && argv[i+2][0]!='-'){
							cout << "ERROR: Only one input file is accepted.\n\n";
							output_help();
							exit(-1);
						}
					}
					break;
				// checking for file with the original filenames, that comes from parameters ini file (cf...)
				case 'f' :
					if(strlen(argv[i]) == 2 || (strlen(argv[i]) > 2 && string(argv[i])=="-file" ) ){
						if((i+1) < argc && argv[i+1][0]!='-'){
							file_filenames = argv[i+1];
						}
						if((i+2) < argc && argv[i+2][0]!='-'){
							cout << "ERROR: Only one filenames file is accepted.\n\n";
							output_help();
							exit(-1);
						}
					}
					break;
				// checking for file with the original filenames, that comes from parameters ini file (cf...)
				case 'p' :
					if(string(argv[i])=="-pre"){
						if((i+1) < argc && argv[i+1][0]!='-'){
							pre_path = argv[i+1];
						}
						if((i+2) < argc && argv[i+2][0]!='-'){
							cout << "ERROR: Only one pre path name is accepted.\n\n";
							output_help();
							exit(-1);
						}
					}
					if(string(argv[i])=="-post"){
						if((i+1) < argc && argv[i+1][0]!='-'){
							post_path = argv[i+1];
						}
						if((i+2) < argc && argv[i+2][0]!='-'){
							cout << "ERROR: Only one pre path name is accepted.\n\n";
							output_help();
							exit(-1);
						}
					}
					break;
				// optional output filename
				case 'o' :
					if(strlen(argv[i]) == 2 || (strlen(argv[i]) > 2 && string(argv[i])=="-out" ) ){
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
	
	// checking if you have input files
	if(file_in.length() == 0){
		cout << "ERROR: Input file is mandatory.\n\n";
		output_help();
		exit(-1);
	}
	
	if(file_filenames.length() == 0){
		cout << "ERROR: Filenames file is mandatory.\n\n";
		output_help();
		exit(-1);
	}
	
	lookbackCat(file_in, file_filenames, pre_path, post_path, file_out);
	
	return 0;
 }