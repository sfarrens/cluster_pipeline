/*
	Code:			convert_hdf5mock.cpp
	Author:         Walter A. Santos Jr. (walter.augusto@gmail.com)
	Last update:	29/09/2012
	Version: 		1.0
	Description: 	Code to convert a hdf5 mock file (from Alex/Euclid) into an ASCII file.
					Makes use of hdf5 lib (http://www.hdfgroup.org/)
	Usage:			see output_help() below
					
	Compile (in VS): cl /MD /EHsc convert_hdf5mock.cpp hdf5mock.cpp /I"C:\Program Files\HDF Group\HDF5\1.8.9\include" /I"C:\Program Files\HDF Group\HDF5\1.8.9\include\cpp" hdf_lib\hdf5_cppdll.lib hdf_lib\hdf5dll.lib
*/

#include <sstream>
#include "hdf5mock.h"

using namespace std;

const char* DEFAULT_VARS[] = {"GalaxyID", "DHaloID", "idrep", "ra", "dec", "z_cos"};


void output_help(){
	cout<<"\nUsage:\n"<<flush;
	cout<<"./hdf5mock -inp <file_in> (-out <file_out>) (-var <column_1>...<column_n>)\n"<<flush;
	cout<<"where <file_in> should be the input HDF5 mock catalog\n"<<flush;
	cout<<"It outputs a set of default variables/columns (id, ra, dec, z...) or given by -var optionally.\n\n"<<flush;
	cout<<"Outputs to <file_out> or to stdout.\n\n"<<flush;
}


int main(int argc, char *argv[]){

	// Checking command line
		if(argc == 1){
		output_help();
		return 0;
	}
	
	// First definitions, before further checking command line
	string file_in, file_out;
	vector<string> vars;
	vector<string> vars_default(DEFAULT_VARS, DEFAULT_VARS + sizeof(DEFAULT_VARS)/sizeof(DEFAULT_VARS[0]));
			
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
				// input hdf5 file check, mandatory option
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
				// output filename, optional
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
				// variables (optional)
				case 'v' :
					if((i+1) < argc && strlen(argv[i]) == 2 || (strlen(argv[i]) > 2 && string(argv[i])=="-var" ) ){
						for (int j=i+1; j<argc; j++){
							if(argv[j][0]=='-')
								break;
							else
								vars.push_back(argv[j]);
						}
					}
					break;
			}
			
		}
	}
	
	// checking for hdf5 file, mandatory
	if(file_in.length() == 0){
		cout << "ERROR: Input hdf5 file is mandatory.\n\n";
		output_help();
		exit(-1);
	}
	
	// checking if vars were given in the input; if not, use default
	if(vars.empty())
		vars = vars_default;
	
	int nvars = vars.size();
	
	vector <vector <string> > data_all;
	
	for(int i=0; i<nvars; i++)
		data_all.push_back(vector<string>());

	string type;
	//stringstream out_temp (stringstream::in | stringstream::out);
	int nrows;
	
	
		
	for(int i=0; i<nvars; i++){
		type = hdf5mock::getType(file_in, "/Data/"+vars[i]);
		
		// gettype types: "int", "long long", "float", "double"
		if(type.compare("int") == 0){
			
			vector<int> data_out = hdf5mock::readColumnInt(file_in, "/Data/"+vars[i]);
			nrows = data_out.size();
			if (nrows){
				for(int j=0; j<nrows; j++){
					stringstream out_temp;
					out_temp << data_out[j];
					data_all[i].push_back(out_temp.str());
					
				}
			}
		}
		else if(type.compare("long long") == 0){
			vector<long long> data_out = hdf5mock::readColumnLongLong(file_in, "/Data/"+vars[i]);
			
			nrows = data_out.size();
			
			
			if (nrows){
				for(int j=0; j<nrows; j++){
					stringstream out_temp;
					out_temp << data_out[j];
					data_all[i].push_back(out_temp.str());
					//if (j==nrows-1)
						//cout << data_all[i][j] << endl;
					
				}
			}
			
		}
		else if(type.compare("float") == 0){
			vector<float> data_out = hdf5mock::readColumnFloat(file_in, "/Data/"+vars[i]);
			nrows = data_out.size();
			if (nrows){
				for(int j=0; j<nrows; j++){
					stringstream out_temp;
					out_temp << data_out[j];
					data_all[i].push_back(out_temp.str());
					
				}
			}
		}
		else if(type.compare("double") == 0){
			vector<double> data_out = hdf5mock::readColumnDouble(file_in, "/Data/"+vars[i]);
			nrows = data_out.size();
			if (nrows){
				for(int j=0; j<nrows; j++){
					stringstream out_temp;
					out_temp << data_out[j];
					data_all[i].push_back(out_temp.str());
					
				}
			}
		}
		
	}
		
	ofstream out;
	if(file_out.length() > 0)
		out.open(file_out.c_str());
		
	//cout << nrows << endl;
	
	for (int i=0; i<nrows; i++){
		for (int j=0; j<nvars; j++){
			if(file_out.length() > 0)
				out << data_all[j][i] << " ";
			else
				cout << data_all[j][i] << " ";
		}
		if(file_out.length() > 0)
			out << "\n";
		else
			cout << "\n";
	}	

	return 0;
}
