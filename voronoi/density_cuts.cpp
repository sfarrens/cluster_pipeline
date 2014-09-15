/*
	Code:			density_cuts.cpp
	Author:         Walter A. Santos Jr. (walter.augusto@gmail.com)
	Last update:	08/08/2013
	Version:		1.1
	Description: 	Code to cut a distribution of galaxies based on voronoi density cuts
					as a function of redsift.
	Last Changes:
					- Bug fixes
					- Makes use of wutils again
					- Fixed output_help()
					- Includes several comments
	Usage:	see output_help() below.
*/

// Libraries
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <sstream> 

#include "wutils.h"

using namespace std;

const double DEFAULT_PERC_CUT = 0.1; // default percentage for density cut, based on redshift bins

//Help / Usage:
void output_help(){
	cout<<"\nUsage:\n"<<flush;
	cout<<"./density_cuts -i <radec_voro_file> (-cut <perc_cut>) (-o <file_out>)\n\n"<<flush;
	cout<<"<radec_voro_file> is the input catalogs\n"<<flush;
	cout<<"with the following columns: id ra dec redshift density (...)\n\n"<<flush;
	cout<<"<perc_cut> is the cutting percentage for each redshift bin (optional, default="<< DEFAULT_PERC_CUT <<".\n"<<flush;
	cout<<"For no cutting, then use '-c 1'.\n\n"<<flush;
	cout<<"Output is a filename, optional, for the output (with same columns as input).\n\n"<<flush;
}

// function to sort points int ids based on sorting redshifts
void sort_with_ids(std::vector<double> &arr, std::vector<int> &brr)
{
	const int M=7,NSTACK=50;
	int i,ir,j,k,jstack=-1,l=0;
	double a;
	int b;
	std::vector<int> istack(NSTACK);

	int n=arr.size();
	ir=n-1;
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j];
				b=brr[j];
				for (i=j-1;i>=l;i--) {
					if (arr[i] <= a) break;
					arr[i+1]=arr[i];
					brr[i+1]=brr[i];
				}
				arr[i+1]=a;
				brr[i+1]=b;
			}
			if (jstack < 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			std::swap(arr[k],arr[l+1]);
			std::swap(brr[k],brr[l+1]);
			if (arr[l] > arr[ir]) {
				std::swap(arr[l],arr[ir]);
				std::swap(brr[l],brr[ir]);
			}
			if (arr[l+1] > arr[ir]) {
				std::swap(arr[l+1],arr[ir]);
				std::swap(brr[l+1],brr[ir]);
			}
			if (arr[l] > arr[l+1]) {
				std::swap(arr[l],arr[l+1]);
				std::swap(brr[l],brr[l+1]);
			}
			i=l+1;
			j=ir;
			a=arr[l+1];
			b=brr[l+1];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				std::swap(arr[i],arr[j]);
				std::swap(brr[i],brr[j]);
			}
			arr[l+1]=arr[j];
			arr[j]=a;
			brr[l+1]=brr[j];
			brr[j]=b;
			jstack += 2;
			if (jstack >= NSTACK){ std::cout << "NSTACK too small in sort2." << std::endl; exit(-1);}
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
}

// function to create redshift histogram bins, and populate them
vector<int> create_histogram(vector<double> redshift_vec, vector<int>& hist_pos){
	int nrows = redshift_vec.size();
	double val_max, val_min;

	wutils::minmax(redshift_vec, nrows, val_min, val_max);
	//minmax(redshift_vec, nrows, val_min, val_max);

//	int num_bins = (int)(val_max - val_min + 1);
	int num_bins = (int) sqrt(nrows);
	vector<int> bins;
	bins.resize(num_bins, 0);
	double bin_width = (val_max-val_min)/num_bins;
	int bin_idx;



	for(int i=0; i<nrows; i++){
		hist_pos[i] = (int)((redshift_vec[i] - val_min - 0.00001) / bin_width);
		if (hist_pos[i] < 0)
			hist_pos[i] = 0;
		bins[hist_pos[i]]++;
	}
	return bins;
}

int main(int argc, char *argv[]){

	// Checking command line
		if(argc == 1){
		output_help();
		return 0;
	}
	
	// First definitions, before further checking command line
	string file_in, file_out;
	double perc_cut = DEFAULT_PERC_CUT;
	
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
				// input voro file (after lookback with ra/dec/redshift) check, mandatory option
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
				// percentage density cut input check, optional;
				case 'c' :
					if(strlen(argv[i]) == 2 || (strlen(argv[i]) > 2 && string(argv[i])=="-cut" )){
						if((i+1) < argc && argv[i+1][0]!='-' && atof(argv[i+1]) >= 0 && atof(argv[i+1]) <= 1){
							perc_cut = atof(argv[i+1]);
						}
					}
					break;
				//output filename, optional
				case 'o' :
					if((i+1) < argc && strlen(argv[i]) == 2 || (strlen(argv[i]) > 2 && string(argv[i])=="-out" ) ){
						if((i+1) < argc && argv[i+1][0]!='-')
							file_out = argv[i+1];
						if((i+2) < argc && argv[i+2][0]!='-'){
							cout << "ERROR: Only one output file is accpeted.\n\n";
							output_help();
							exit(-1);
						}
					}
					break;
			}
			
		}
	}
		
	// checking for input file, mandatory
	if(file_in.empty()){
		cout << "ERROR: Input file(s) is mandatory.\n\n";
		output_help();
		exit(-1);
	}
	// assuming...
	// 0,1,2,3,4
	//id, ra, dec, redshift, density
	int ncols, nrows;
	vector< vector<string> > data_mx;
	wutils::read_ascii(file_in, data_mx, nrows, ncols);
	//read_ascii(file_in, data_mx, nrows, ncols);
	
	vector<double> redshifts;
	redshifts.resize(nrows);
	vector<double> densities;
	densities.resize(nrows);

	//assuming redshift column #3, density #4
	for(int i=0; i<nrows; i++){
		redshifts[i] = atof(data_mx[3][i].c_str());
		densities[i] = atof(data_mx[4][i].c_str());
	}
	vector<int> bin_counts;
	vector<int> hist_positions;
	hist_positions.resize(nrows);

	// create histogram from redshifts and also gives bin position for each point
	bin_counts = create_histogram(redshifts, hist_positions);
	int nbins = bin_counts.size();
	vector< vector <int> > p_ids;
	p_ids.resize(nbins);
	vector< vector <double> > p_densities;
	p_densities.resize(nbins);

	// "populating" bins with ids and densities
	for(int i=0; i<nrows; i++){
		p_ids[hist_positions[i]].push_back(i);
		p_densities[hist_positions[i]].push_back(densities[i]);
	}
	vector<int> surviving_ids;
	int npoints = 0;
	int n_surviving = 0;

	// sorting desnities in each bin and applying the cut
	for(int k=0; k<nbins; k++){
		npoints = bin_counts[k];
		vector<double> k_densities;
		vector<int> k_ids;
		for(int f=0; f<npoints; f++){
			k_densities.push_back(p_densities[k][f]);
			k_ids.push_back(p_ids[k][f]);
		}

		sort_with_ids(k_densities, k_ids);

		n_surviving = (int)floor(perc_cut*npoints);
		for(int j=0; j<n_surviving; j++)
			surviving_ids.push_back(k_ids[j]);
	}

	ofstream out;
	if(file_out.length() > 0)
		out.open(file_out.c_str());
			
	//coping the surviving data
	int n_surviving_total = surviving_ids.size();

	for(int j=0; j<n_surviving_total; j++){
		for(int c=0; c<ncols; c++){
			if(file_out.length() > 0)
				out << data_mx[c][surviving_ids[j]] << "  ";
			else
				cout << data_mx[c][surviving_ids[j]] << "  ";
		}
		if(file_out.length() > 0)
			out << "\n";
		else
			cout << "\n";
	}
	
	return 0;
}
