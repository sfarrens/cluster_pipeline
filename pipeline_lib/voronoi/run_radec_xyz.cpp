/*
	Code:			run_radec_xyz.cpp
	Author:         Walter A. Santos Jr. (walter.augusto@gmail.com)
	Last update:	14/08/2012
	Version: 		1.2
	Last changes:	
		- Fixed a bug with the out filename;
		- Included cosm option directly in the command line; 
		- cout only outputs the results now, unless it's an error/exit
	Description: 	Code to transform from radecz into xyz space, from a catalog of objects.
					Implements class radec_xyz in radec_xyz.h. Makes use of wutils.h.
	Compile (in VS!): cl /MD /EHsc run_radec_xyz.cpp cfit3300\cfitsio.lib radec_xyz.cpp wutils.cpp

	Usage:			see output_help() function below.
*/

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

#include "radec_xyz.h"
#include "wutils.h"

using namespace std;

const double DEFAULT_OMEGA_M = 0.3;
const double DEFAULT_HUBBLE = 0.7;
const double DEFAULT_OMEGA_L = 0.7;

void output_help(){
	cout<<"\nUsage:\n"<<flush;
	cout<<"./run_radec_xyz -inp <file_in> (-out <file_out>) (-om <%f>) (-ol <%f>) (-hubble <%f>)\n"<<flush;
	cout<<"(-distfile <file_dist>) (-param <param.cosm>)\n\n"<<flush;
	cout<<"where <file_in> should be the input catalog\n"<<flush;
	cout<<"with the following columns: id ra dec redshift (...)\n\n"<<flush;
	cout<<"It outputs 'id x y z' to stdout or to file <file_out>.\n\n"<<flush;
	cout<<"There are 3 (non-mandatory) option to include the cosmology (omega_l, omega_m, h):\n"<<flush;
	cout<<"1) Directly to the command line, with options: -om -ol -hubble. It overrides values in option 2.\n"<<flush;
	cout<<"2) Param file with cosmological parameters: omega_l, omega_m, h: -param.\n."<<flush;
	cout<<"3) Input a file with distances as a function of 'a' (a x distance) : -distfile. It overrides options 1 and 2.\n"<<flush;
	cout<<"If none of these are provided (or if the values are wrong), following defaults are used:\n"<<flush;
	cout<<"omega_l = " << DEFAULT_OMEGA_L << " ; omega_m = " << DEFAULT_OMEGA_M << " ; hubble = " << DEFAULT_HUBBLE << "\n\n"<<flush;
}

void read_cosm(string file_param, double& omega_m, double& omega_l, double& hubble_h){
	vector < vector <string> > params;
	int nrows, ncols;
	size_t found_om, found_ol, found_h;
	
	wutils::read_ascii(file_param, params, nrows, ncols);
	
	for(int i=0; i<nrows; i++){
		found_om = params[0][i].find("omega_m");
		found_ol = params[0][i].find("omega_l");
		found_h = params[0][i].find("h");
		
		if (found_om!=string::npos){
			//cout << "found_om" << endl;
			omega_m = atof(params[1][i].c_str());
		}
		else if (found_ol!=string::npos){
			//cout << "found_ol" << endl;
			omega_l = atof(params[1][i].c_str());
		}
		else if (found_h!=string::npos){
			//cout << "found_h" << endl;
			hubble_h = atof(params[1][i].c_str());
		}
	}
	// cout << "omega_m: " << omega_m << endl;
	// cout << "omega_l: " << omega_l << endl;
	// cout << "h: " << hubble_h << endl;
}

int main(int argc, char *argv[]){
	string file_in, file_out, file_dist, file_param;
	double omega_m=-1, omega_l=-1, hubble=-1;

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
					// checking for hubble parameter option
					else if((i+1) < argc && strlen(argv[i]) > 2 && string(argv[i])=="-hubble"){
						hubble = atof(argv[i+1]);
					}
					break;
				// checking for input file	
				case 'i' :
					if((i+1) < argc && (strlen(argv[i]) == 2 || (strlen(argv[i]) > 2 && string(argv[i])=="-inp" )) ){
						file_in = argv[i+1];
						if((i+2) < argc && argv[i+2][0]!='-'){
							cout << "ERROR: Only one input file is accepted.\n\n";
							output_help();
							exit(-1);
						}
					}
					break;
				
				case 'o' :
					if((i+1) < argc){
						// checking for output file
						if(strlen(argv[i]) == 2 || (strlen(argv[i]) > 2 && string(argv[i])=="-out" ) ){
							file_out = argv[i+1];
							if((i+2) < argc && argv[i+2][0]!='-'){
								cout << "ERROR: Only one output file is accepted.\n\n";
								output_help();
								exit(-1);
							}
						}
						// checking for omega_m
						else if(string(argv[i])=="-om" || string(argv[i])=="-omega_m"){
							omega_m = atof(argv[i+1]);
						}
						// checking for omega_l
						else if(string(argv[i])=="-ol" || string(argv[i])=="-omega_l"){
							omega_l = atof(argv[i+1]);
						}
					}
					break;
				//checking for cosm distance x a file
				case 'd' :
					if((i+1) < argc && string(argv[i])=="-distfile")
						file_dist = argv[i+1];
					break;
				// checking for param cosm file
				case 'p' :
					if((i+1) < argc && string(argv[i])=="-param")
						file_param = argv[i+1];
					break;
			}
			
		}
	}
	
	// checking if you have input file
	if(file_in.empty()){
		cout << "ERROR: Input file is mandatory.\n\n";
		output_help();
		exit(-1);
	}
	
	
	
	// use of the class radec_xyz (defined in radec_xyz.h) depends on the input...
	radec_xyz coord;
	
	if(!(file_dist.empty())) // reading distances x a from a file, overrides anything else
		coord.setup(file_dist);
	else {
		//go through each cosm parameter
		double om, ol, h; //aux params
		
		if (!(file_param.empty())) // try to reading cosmology (omega_m, omega_l, h) from a param file
				read_cosm(file_param, om, ol, h);
		
		//checking omega_m
		if(omega_m < 0 || omega_m > 1){ // it means the input -om not used or is wrong value
			if (!(file_param.empty())){ // from the param file
				if (om >= 0 && om <=1)
					omega_m = om;
				else
					omega_m = DEFAULT_OMEGA_M;
			}
			else
				omega_m = DEFAULT_OMEGA_M;			
		}
		
		
		
		//checking omega_l
		if(omega_l < 0 || omega_l > 1){ // it means the input -ol not used or is wrong value
			if (!(file_param.empty())){ // from the param file
				if (ol >= 0 && ol <=1)
					omega_l = ol;
				else
					omega_l = DEFAULT_OMEGA_L;
			}
			else
				omega_l = DEFAULT_OMEGA_L;			
		}
	
		//checking h hubble
		if(hubble < 0 || hubble > 1){ // it means the input -h not used or is wrong value
			if (!(file_param.empty())){ // from the param file
				if (h >= 0 && h <=1)
					hubble = h;
				else
					hubble = DEFAULT_HUBBLE;
			}
			else
				hubble = DEFAULT_HUBBLE;			
		}
	
		
		// setting up with the values
		coord.setup(omega_m, omega_l, hubble);
	}
			
	double x, y, z;
	int ncols, nrows;
	vector< vector<string> > data_mx;
	
	//wutils::read_fits(file_in, data_mx, nrows, ncols);
		
	wutils::readFile(file_in, data_mx, nrows, ncols);
	
	
	ofstream out;
	if(file_out.length() > 0)
	// file_out[k] = wutils::tag_outfile(file_in[k], tag);
		out.open(file_out.c_str());
		
	for (int i=0; i < nrows; i++){
		coord.xyz(atof(data_mx[1][i].c_str()), atof(data_mx[2][i].c_str()), atof(data_mx[3][i].c_str()), x, y, z);
		if(file_out.length() > 0)
		  out << data_mx[0][i] << " " << x << " " << y << " " << z << " " << file_in << endl;
		else
		  cout << data_mx[0][i] << " " << x << " " << y << " " << z << " " << file_in << endl;
	}
	if(file_out.length() > 0){
		cout << "Output saved in " << file_out << "." << endl;
		out.close();
	}
	
	return 0;
}
