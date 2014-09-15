/*
	Code:			main_split_radec.cpp
	Author:         Walter A. Santos Jr. (walter.augusto@gmail.com)
	Last update:	15/08/2013
	Version:		1.1
	Description: 	Code to split a radec space into boxes, from a list of object positions/files.
					It must be run *before* run_split_radec. Requires the 'bound' file.
	Last Changes:	- It accepts number of boxes instead of average number of points in each box.
					(It doesn't read anymore all input files for the total number of points)
	Compile:		
	Usage:	see output_help() below.
*/

// Libraries
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

#include "wutils.h"

using namespace std;

const int NUM_DIMS=2; // radec
//const int DEFAULT_NUM_MAX_CUT=5000; // default average number of data points in each box/file
const int DEFAULT_NUM_BOXES=4; // default average number of data points in each box/file
const double DEFAULT_INTERS_PERC = 0.05; // default percentage of intersection between boxes

// structure for a box
struct box {
	double split_v; // splitting value in the split_d dimension of the node
	short split_d; // splitting dimension of the node
	int rank;	// how far away from the root/overall box
	short leaf;  // flag for leaf, 1 means is a leaf and will not be further split
	std::vector<double> min, max; // min/max in each dimension
	struct box *left, *right; // recursive child boxes
};

// save all leaves (and only leaves) recursively into vectors
void save_leaf_boxes(box *node, vector< vector<double> > &min, vector< vector<double> > &max){
	if (node->leaf){
		for(int j=0; j<NUM_DIMS; j++){
			min[j].push_back(node->min[j]);
			max[j].push_back(node->max[j]);
		}
	}
    else {
        save_leaf_boxes(node->left, min, max);
        save_leaf_boxes(node->right, min, max);
    }
}

// Function to find the widest dimension of a box (splitting dimension)
int widest_dim(std::vector<double> min, std::vector<double> max, int dims){
    int wd=1, k;
    double d, dmax=0;
        
    for (k=0; k<dims; k++){
        d = max[k] - min[k];
        if (d > dmax){
            dmax = d;
            wd = k;
        }
    }
    return wd;
}

// Function to create the boxes recursively
box *create_boxes(const int num_dims, const int max_split_times, int rank, std::vector<double> min, std::vector<double> max,
					double inters){
	box *node = new box;
	
	node->min = min;
	node->max = max;
	node->rank = rank;
	
	if (rank >= max_split_times) { // leaf test
	    node->leaf = 1;
		node->left = NULL;
		node->right = NULL;
		node->split_d = -1;
		node->split_v = 0;
        return node;        
    }
	
	// if not leaf, keep on splitting
	node->leaf = 0;
	
	node->split_d = widest_dim(min, max, num_dims); // splitting dimension is the widest
	node->split_v = (max[node->split_d]+min[node->split_d])/2; // splitting value is the mean point
	double width = max[node->split_d]-min[node->split_d]; // width in the splitting dimension
	
	std::vector<double> min_left = min;
	std::vector<double> max_left = max;
	// setting intersection space, as a percentage of the width in the splitting dimension
	max_left[node->split_d] = node->split_v + width*inters;
	std::vector<double> min_right = min;
	min_right[node->split_d] = node->split_v - width*inters; //again, setting intersection space
	std::vector<double> max_right = max;
	
	// recursively create the sub-boxes
	node->left = create_boxes(num_dims, max_split_times, rank+1, min_left, max_left, inters);
	node->right = create_boxes(num_dims, max_split_times, rank+1, min_right, max_right, inters);
	    
    return node;
}

void output_help(){
	cout<<"\nUsage:\n"<<flush;
	cout<<"./main_split -b <file_boundary> (-inters <perc>) (-n <N>) (-o <file_out>)\n\n"<<flush;
	cout<<"<file_boundary> is the file containing the overall boundaries (min/max)\n"<<flush;
	cout<<"<perc> is the percentage of intersection between boxes (optional, default="<<DEFAULT_INTERS_PERC<<")\n"<<flush;
	cout<<"<N> is the number of boxes to divide into (optional, default="<<DEFAULT_NUM_BOXES<<")\n\n"<<flush;
	cout<<"<file_out> is the output filename; if not given, default it to stdout.\n\n"<<flush;
}


// main

int main(int argc, char *argv[]){

	// Checking command line
		if(argc == 1){
		output_help();
		return 0;
	}
	
	// First definitions, before further checking command line
	vector<string> file_in;
	string file_bound, file_boxes;
	double inters = DEFAULT_INTERS_PERC;
	//int num_max_split = DEFAULT_NUM_MAX_CUT;
	int num_boxes_input = DEFAULT_NUM_BOXES;
		
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
						file_bound = argv[i+1];
						if((i+2) < argc && argv[i+2][0]!='-'){
							cout << "ERROR: Only one input file is accepted.\n\n";
							output_help();
							exit(-1);
						}
					}
					break;
				// input xyz file check, mandatory option
				/*
				case 'i' :
					if((i+1) < argc){
						if(strlen(argv[i]) == 2 || (strlen(argv[i]) > 2 && string(argv[i])=="-inp" )){
							for (int j=i+1; j<argc; j++){
								if(argv[j][0]=='-')
									break;
								else
									file_in.push_back(argv[j]);
							}
						}
						// optional intersection parameter
						else if(strlen(argv[i]) > 2 && string(argv[i])=="-inters"){
							if(argv[i+1][0]!='-' && atof(argv[i+1]) >= 0 && atof(argv[i+1]) < 1){
								inters = atof(argv[i+1]);
								
							}
						}
					}
					break;
					*/
				// average data points input check, optional; default = ?
				case 'n' :
					if(strlen(argv[i]) == 2){
						if((i+1) < argc && argv[i+1][0]!='-' && atoi(argv[i+1]) > 0){
							num_boxes_input = atoi(argv[i+1]);
						}
					}
					break;
				case 'o' :
					if((i+1) < argc && (strlen(argv[i]) == 2 || (strlen(argv[i]) > 2 && string(argv[i])=="-out"))){
						// checking for output file
						file_boxes = argv[i+1];
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
	/*
	// checking for xyz file, mandatory
	if(file_in.empty()){
		cout << "ERROR: Input file(s) is mandatory.\n\n";
		output_help();
		exit(-1);
	}
	*/
	
	//int nfiles = file_in.size();
	//cout << "Number of input xyz files: " << nfiles << endl;
	//cout << "Intersection percentage: " << inters << endl;
	
	// reading and saving boundary data
	std::vector<double> min_all;
	std::vector<double> max_all;
	min_all.resize(NUM_DIMS);
	max_all.resize(NUM_DIMS);
	int ncols, nrows;
	vector< vector<string> > data_mx;
	
	wutils::read_ascii(file_bound, data_mx, nrows, ncols);
	for(int k=0; k<NUM_DIMS; k++){
		min_all[k] = atof(data_mx[0][k].c_str());
		max_all[k] = atof(data_mx[1][k].c_str());
	}
	/*
	//counting all data points
	long ndata=0;
	for (int k=0; k<nfiles; k++){
		vector< vector<string> > count_mx;
		wutils::read_ascii(file_in[k], count_mx, nrows, ncols);
		ndata += nrows;
	}
	
	// check if spplitting is needed, given the average/max data points
	if (ndata <= num_max_split){
		cout << "No splitting is needed. Input data points less than " << num_max_split <<"."<<endl;
		cout << "Proceed to step voro++." << endl;
		return 0;
	}
	*/
	
	// define the number of boxes (in double for now) and how many times of the splliting
	//double num_boxes = (double) ndata / (double) num_max_split;
	//int max_splitting = (int) ceil(log(num_boxes)/log(2.0));
	int max_splitting = (int) ceil(log(num_boxes_input)/log(2.0));
	
	// create all boxes, dividing the space
	box *allboxes = create_boxes(NUM_DIMS, max_splitting, 0, min_all, max_all, inters);
	
	// save the sub-boxes into min/max vectors
	vector< vector<double> > min_boxes;
	vector< vector<double> > max_boxes;
	for(int j=0;j<NUM_DIMS;j++){
		min_boxes.push_back(vector<double>());
		max_boxes.push_back(vector<double>());
	}
	save_leaf_boxes(allboxes, min_boxes, max_boxes);
	
	// official number of boxes/leaves
	//num_boxes = min_boxes[0].size();
	int num_boxes = min_boxes[0].size();
	//cout << "Number of splitting boxes: " << num_boxes << endl;
	//cout << "\n";
	
	// saving a file with the out filenames, min, max for each sub-box
	//string file_boxes = wutils::replace_tag(file_in[0], "xyz", "boxes");
	ofstream out;
	if(file_boxes.length() > 0)
		out.open(file_boxes.c_str());
	
	for (int i=0; i<num_boxes; i++){
		//file_out = wutils::tag_outfile(file_in[0], "split."+wutils::IntToStr(i));
		if(file_boxes.length() > 0){
			out << i << " ";
			for (int j=0; j<NUM_DIMS; j++)
				out << min_boxes[j][i] << " " << max_boxes[j][i] << " ";
			out << "\n";
		}
		else{
			cout << i << " ";
			for (int j=0; j<NUM_DIMS; j++)
				cout << min_boxes[j][i] << " " << max_boxes[j][i] << " ";
			cout << "\n";
		}
	}
	
	if(file_boxes.length() > 0){
		cout << "Output saved in " << file_boxes << "." << endl;
		out.close();
	}
		
	return 0;
}
