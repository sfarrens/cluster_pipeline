/*
	Code:			radec_xyz.h
	Author:         Walter A. Santos Jr. (walter.augusto@gmail.com)
	Last update:	12/07/2012
	Description: 	Converts ra,dec,redshift into x,y,z coordinates.
	Usage:			Header/Library file. Include "radec_xyz.h" to use in other C/C++ codes.
					Source codes are in radec_xyz.cpp
					radec should be given in degrees.
					Example:
						radec_xyz coord;	// default constructor -> default cosmology constants
						double ra = 100, dec = 10, redshift = 0.2;		// radec in degress
						double x, y, z;
						coord.xyz(ra, dec, redshift, x, y, z);		// xyz are given the converted values
*/

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <cfloat>
#include <cstdlib>
#include <cstring>

#include "wutils.h"

#ifndef RADEC_XYZ_H
#define RADEC_XYZ_H

class radec_xyz
{
	public:
		// constructor
		// uses default cosm constants if none are given
		radec_xyz(double omega_m=0.3, double omega_l=0.7, double hubble_h=0.7) : om(omega_m),
					ol(omega_l), h(hubble_h), iscalcdm(1) { }
		radec_xyz(std::string distfile) : _distfile(distfile), iscalcdm(0) { setup_infile(); }
		// if the default constructor was already used, this will setup the class again
		// with dist file, or new values for cosm, etc.
		void setup (double omega_m, double omega_l, double hubble_h){
			om = omega_m;
			ol = omega_l;
			h = hubble_h;
			iscalcdm = 1;
		}
		void setup (std::string distfile){
			_distfile = distfile;
			iscalcdm = 0;
			setup_infile();
		}
		// functions to calculate xyz from ra,dec,redshift; radec in degress
		// calculates and returns only x
		double x(double ra, double dec, double redshift);
		// calculates and returns only y
		double y(double ra, double dec, double redshift);
		// calculates and returns only z
		double z(double ra, double dec, double redshift);
		// calculates xyz, updating given xyz variables
		void xyz(double ra, double dec, double redshift, double& _x, double& _y, double& _z);
	private:
		std::string _distfile;
		std::vector < std::vector < std::string > > data_mx;
		std::vector <double> inp_distances, inp_redshifts, y2_splines;
		int nrows, ncols;
		short iscalcdm;
		double om, ol, h;
		double f_dist(double z);
		double trapzd(double a, double b, int n);
		double calc_dm(double z);
		double calc_x(double ra, double dec, double dm);
		double calc_y(double ra, double dec, double dm);
		double calc_z(double dec, double dm);
		double redshiftFromA(double a){
			if (a == 0)
				return -1;
			else
				return 1/a - 1;
		}
		void setup_infile();
		double lookup_distance(double redshift);
};

#endif // RADEC_XYZ_H
