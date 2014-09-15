/*
	Code:			utils.h
	Author:         Walter A. Santos Jr. (walter.augusto@gmail.com)
	Based on:		Some functions based on Sam Farrens's and UCLs Weak Lensing Pipeline (WLP)
					via Filipe Abdalla
	Last update:	24/07/2012
	Description: 	Utility codes to use on other routines, like Voronoi3D.
	Usage:			Header/Library file. Include "utils.h" to use in other C/C++ codes.
					Call functions with wutils::function
*/

#ifndef WUTILS_H
#define WUTILS_H

// libraries
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>   
#include <cmath>
#include <algorithm>

//#include "../cfitsio/longnam.h"
//#include "../cfitsio/fitsio.h"

namespace wutils {

		// void SWAP(double a, double b){
			// double temp = a;
			// a = b;
			// b = temp;
		// }		
		void sort2(std::vector<double> &arr, std::vector<double> &brr);
		void spline(std::vector<double> &x, std::vector<double> &y, const double yp1, const double ypn,
	std::vector<double> &y2);
	
		void splint(std::vector<double> &xa, std::vector<double> &ya, std::vector<double> &y2a, const double x, double &y);
		
		std::string IntToStr(int n);
		
		// given an array of double values with size num_points, find min and max values
		void minmax(std::vector<double>& list, const int num_points, double& min, double& max);

		// given a 3d matrix of double values, find min and max boundaries as vectors
		void find_boundaries(std::vector< std::vector< double > > &input_mx, const int num_dims, const int num_points, std::vector<double>& min, std::vector<double>& max);

		// given an input filename string, return an output filename with a suffix to be appended
		std::string tag_outfile(const std::string &in, const std::string &tag);
		
		std::string replace_tag(const std::string &in, const std::string &oldtag, const std::string &newtag);
		
		std::string check_file_type(const std::string &fname);
		
		int count_cols(const std::string &fname);
		
		void skip_comments(std::ifstream& fstr);
		
		void read_ascii(const std::string &fname, std::vector< std::vector< std::string > > &data_mx, int& nrows, int& ncols);
		
		//void read_fits(const std::string &fname, std::vector< std::vector< std::string > > &data_mx, int& nrows, int& ncols);
		
		// read a file, saving the data into a matrix rows x col, nrows, ncols)
		//void readFile(const std::string& filename, std::vector< std::vector< std::string > > &data_mx, int& nrows, int& ncols);

		// alloc/free memory functions, type double
		double **alloc_2dmatrix(int rows, int cols);
		void free_2dmatrix(double **mx, int rows);
		double ***alloc_3dmatrix(int x, int y, int z);
		void free_3dmatrix(double ***mx, int x, int y);
}


#endif // WUTILS_H
