/*
	Code:			find_boundaries_radec.cpp
	Author:         Walter A. Santos Jr. (walter.augusto@gmail.com)
	Last update:	07/08/2013
	Version: 		0.1
	Description: 	Code to to find boundaries in radec space, from a list of object positions.
					Makes use of wutils.h.
	Usage:			see output_help() function below.
*/

// Libraries
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

#include "wutils.h"

using namespace std;

const int NUM_DIMS=2; // ra, dec

void output_help(){
	cout<<"\nUsage:\n"<<flush;
	cout<<"./find_boundaries -inp <file_1> ... (<file_n>) (-out <file_out>)\n\n"<<flush;
	cout<<"where <file_1> ... (<file_n>) should be the input catalogs\n"<<flush;
	cout<<"with the following columns: id x y z (...)\n\n"<<flush;
	cout<<"It runs with at least 1 file (n=1)\n"<<flush;
	cout<<"Output is either to <file_out> or to stdout.\n\n"<<flush;
}

// main

int main(int argc, char *argv[]){

	// Checking command line
	if(argc == 1){
		output_help();
		return 0;
	}
	
	vector<string> file_in;
	string file_out;
		
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
					if((i+1) < argc && (strlen(argv[i]) == 2 || (strlen(argv[i]) > 2 && string(argv[i])=="-inp" )) ){
						for (int j=i+1; j<argc; j++){
							if(argv[j][0]=='-')
								break;
							else
								file_in.push_back(argv[j]);
						}
					}
					break;
				case 'o' :
					if((i+1) < argc && strlen(argv[i]) == 2 || (strlen(argv[i]) > 2 && string(argv[i])=="-out" ) )
						file_out = argv[i+1];
					break;
			}
			
		}
	}
	
	if(file_in.empty()){
		cout << "ERROR: Input file(s) is mandatory.\n\n";
		output_help();
		exit(-1);
	}
	int nfiles = file_in.size();
	//cout << "nfiles: " << nfiles << endl;
		
	int ncols, nrows;
	std::vector<double> min_all;
	std::vector<double> max_all;
	min_all.resize(NUM_DIMS);
	max_all.resize(NUM_DIMS);
	
	for (int k=0; k<nfiles; k++){
		vector< vector<string> > data_mx;
		
		//wutils::readFile(file_in[k], data_mx, nrows, ncols);
		wutils::read_ascii(file_in[k], data_mx, nrows, ncols);
		
		vector< vector< double > > pos_mx;
		for(int j=0;j<NUM_DIMS;j++){
			pos_mx.push_back(vector<double>());
			pos_mx[j].resize(nrows);
		}
		for (int i=0; i<nrows; i++){
			for (int j=0; j<NUM_DIMS; j++)
				pos_mx[j][i] = atof(data_mx[j+1][i].c_str());
		}
		
		std::vector<double> min;
		std::vector<double> max;
		min.resize(NUM_DIMS);
		max.resize(NUM_DIMS);
		
		wutils::find_boundaries(pos_mx, NUM_DIMS, nrows, min, max);
		
		// refresh overall min/max (for all files)
		
		if(k==0){
			min_all = min;
			max_all = max;
		}
		else {
			for (int j=0; j<NUM_DIMS; j++){
				if (min[j] < min_all[j])
					min_all[j] = min[j];
				if (max[j] > max_all[j])
					max_all[j] = max[j];
			}
		}
	}
	
	ofstream out;
	if(file_out.length() > 0)
		out.open(file_out.c_str());
		
	for (int k=0; k<NUM_DIMS; k++){
		if(file_out.length() > 0)
			out << min_all[k]-0.1 << " " << max_all[k]+0.1 << endl;
		else
			cout << min_all[k]-0.1 << " " << max_all[k]+0.1 << endl;
	}
  
	//cout << "\n";
	 
	if(file_out.length() > 0){
		out.close();
		cout << "Output saved in " << file_out << "." << endl;
	}
	return 0;
}