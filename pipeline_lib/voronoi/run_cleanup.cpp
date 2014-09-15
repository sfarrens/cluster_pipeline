/*
	Code:			run_cleanup.cpp
	Author:         Walter A. Santos Jr. (walter.augusto@gmail.com)
	Last update:	16/08/2012
	Version: 		1.0
	Description: 	Code to cleanup duplicates (intersecting region) in voronoi results files. Makes use of wutils.h.
	
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
	cout<<"./run_cleanup -inp <file1> <file2>\n\n"<<flush;
	cout<<"where <file1> and <file2> are the input catalogs\n"<<flush;
	cout<<"with the following columns: id x y z orig_file orig_line density num_neighbors ids_neighbors(...)\n\n"<<flush;
	cout<<"No new output is generated: the input files are updated with the cleanup duplicated lines (if found).\n\n"<<flush;
}

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

void cleanCats(std::string file1, std::string file2){

	int nrows1 = count_num_points(file1);
	int nrows2 = count_num_points(file2);
	
	std::ifstream f1 (file1.c_str());
	std::ifstream f2 (file2.c_str());

	if (!(f1.is_open()) || !(f2.is_open())){
		std::cout << "ERROR: cannot open files "<< file1 << " and/or " << file2 << "to cleanup!" << std::endl;
		exit(-1);
	}
	
	vector<int> flag1, flag2; // 0 for no-match, 1 for match and use file1, 2 for match and use file2
	flag1.resize(nrows1, 0);
	flag2.resize(nrows2, 0);
	
	std::string id1, id2, aux;
	std::string line;
		int n1, n2;
	int points1, points2;
	int nwalls1, nwalls2, num_same=0, num_final=0;
	
	// cout << "Num file 1 orig: " << nrows1 << endl;
	// cout << "Num file 2 orig: " << nrows2 << endl;
	
	
	for (int i=0; i<nrows1; i++){
		f1 >> id1;
		for (int j=0; j<nrows2; j++){
			f2 >> id2;
			if (id1.compare(id2) == 0){
				num_same++;
				for(int k=0; k<6; k++){ // reading positions, orig_file, orig_line, and volume
					f1 >> aux;
					f2 >> aux;
				}
				f1 >> n1; // reading number of neighbors
				f2 >> n2;
				nwalls1=0;
				nwalls2=0;
				for (int k=0; k<n1; k++){
					f1 >> points1;
					if (points1<0) // negative neighbor ids == 0
						nwalls1++;
				}
				for (int k=0; k<n2; k++){
					f2 >> points2;
					if (points2<0)
						nwalls2++;
				}
				if (nwalls1 <= nwalls2){
					flag1[i]=1;
					flag2[j]=1;
				}
				else {
					flag1[i]=2;
					flag2[j]=2;
				}
				break;
			}
			getline (f2, line);
		}
		f2.seekg (0, std::ios::beg);
		if (flag1[i]==0)
			getline (f1, line);
		
		
	}
	
	//cout << "Num same: " << num_same << endl;
	
	std::string file1_temp = file1 + ".temp";
	std::string file2_temp = file2 + ".temp";

	std::ofstream f1_out (file1_temp.c_str());
	std::ofstream f2_out (file2_temp.c_str());

	if (!(f1_out.is_open()) || !(f2_out.is_open())){
		std::cout << "ERROR: cannot open output temp files "<< file1_temp << " and/or " << file2_temp << std::endl;
		exit(-1);
	}
					
	f1.seekg (0, std::ios::beg);
	for (int i=0; i<nrows1; i++){
		getline (f1, line);
		if (flag1[i]==0 || flag1[i]==1){
			num_final++;
			f1_out << line << std::endl;
		}
	}
	
	//cout << "Num final 1: " << num_final << endl;
		
	num_final = 0;
	f2.seekg (0, std::ios::beg);
	for (int j=0; j<nrows2; j++){
		getline (f2, line);
		if (flag2[j]==0 || flag2[j]==2){
			num_final++;
			f2_out << line << std::endl;
		}
	}
	
	//cout << "Num final 2: " << num_final << endl;
	
	f1.close();
	f2.close();
	f1_out.close();
	f2_out.close();
	
	remove(file1.c_str());
	remove(file2.c_str());
	
	rename(file1_temp.c_str(), file1.c_str());
	rename(file2_temp.c_str(), file2.c_str());

}

int main(int argc, char *argv[]){
	string file_in1, file_in2;
	
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
				// checking for input files	
				case 'i' :
					if((i+2) < argc && (strlen(argv[i]) == 2 || (strlen(argv[i]) > 2 && string(argv[i])=="-inp" )) ){
						if(argv[i+1][0]!='-')
							file_in1 = argv[i+1];
						if(argv[i+2][0]!='-')
							file_in2 = argv[i+2];
						
						if((i+3) < argc && argv[i+3][0]!='-'){
							cout << "ERROR: Only 2 input files are accepted.\n\n";
							output_help();
							exit(-1);
						}
					}
					break;
			}
			
		}
	}
	
	// checking if you have input files
	if(file_in1.length() == 0 || file_in2.length() == 0){
		cout << "ERROR: 2 input files are mandatory.\n\n";
		output_help();
		exit(-1);
	}
	
	
	
	cleanCats(file_in1, file_in2);
	
	return 0;
 }