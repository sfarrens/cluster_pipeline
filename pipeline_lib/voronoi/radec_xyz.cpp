/*
	Code:			radec_xyz.cpp
	Author:         Walter A. Santos Jr. (walter.augusto@gmail.com)
	Last update:	12/07/2012
	Description: 	Converts ra,dec,redshift into x,y,z coordinates.
	Usage:			Include and read "radec_xyz.h" to use in other C/C++ codes.
*/

#include "radec_xyz.h"

using namespace std;

double radec_xyz::f_dist(double z){
	double res = 1/sqrt(om*(1+z)*(1+z)*(1+z)+ol);
	return res;
}

double radec_xyz::trapzd(double a, double b, int n){
	double x, tnm, sum, del;
	static double s;
	int i, j;
	if (n==1){
		return (s=0.5*(b-a)*(f_dist(a)+f_dist(b)));
	}
	else {
		for (i=1,j=1;j<n-1;j++)
			i <<=1;
		tnm=i;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=i;j++,x+=del)
			sum += f_dist(x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}

double radec_xyz::calc_dm(double z){
	if(iscalcdm){
		int i;
		double dc, dh, dm;
		dh = 3000/(h);
		for (i=1; i<=10; i++)
			dc = trapzd (0, z, i);
		dm = dc * dh;
		return dm;
	}
	else
		return lookup_distance(z);
}

double radec_xyz::x(double ra, double dec, double redshift){
	double dm = calc_dm(redshift);
	return calc_x(ra, dec, dm);
}

double radec_xyz::calc_x(double ra, double dec, double dm){
	double x;
	ra = ra*M_PI/180;
	dec = (90-dec)*M_PI/180;
	x = dm * sin(dec)*cos(ra);
	return x;
}

double radec_xyz::y(double ra, double dec, double redshift){
	double dm = calc_dm(redshift);
	return calc_y(ra, dec, dm);
}

double radec_xyz::calc_y(double ra, double dec, double dm){
	double y;
	ra = ra*M_PI/180;
	dec = (90-dec)*M_PI/180;
	y = dm * sin(dec)*sin(ra);
	return y;
}

double radec_xyz::z(double ra, double dec, double redshift){
	double dm = calc_dm(redshift);
	return calc_z(dec, dm);
}

double radec_xyz::calc_z(double dec, double dm){
	double z;
	dec = (90-dec)*M_PI/180;
	z = dm * cos(dec);
	return z;
}
	
void radec_xyz::xyz(double ra, double dec, double redshift, double& _x, double& _y, double& _z){
	double dm;
	dm = calc_dm(redshift);
	_x = calc_x(ra, dec, dm);
	_y = calc_y(ra, dec, dm);
	_z = calc_z(dec, dm);
}

double radec_xyz::lookup_distance(double redshift){
	if (redshift < 0)
		return -1;
		
	double dist;

	wutils::splint(inp_redshifts, inp_distances, y2_splines, redshift, dist);
	return dist;
	
	// if (redshift < 0)
		// return NULL;
	// int match=0;
	// double diff_current, diff=9999;
	// for (int i=0; i<nrows; i++){
		// if (inp_redshifts[i] >= 0){
			// diff_current = abs(redshift - inp_redshifts[i]);
			// if(diff_current < diff){
				// diff = diff_current;
				// match = i;
			// }
		// }
	// }
	// return inp_distances[match];
}

void radec_xyz::setup_infile(){
	wutils::readFile(_distfile, data_mx, nrows, ncols);
	for (int i=0; i<nrows; i++){
		inp_redshifts[i] = redshiftFromA(atof(data_mx[i][0].c_str()));
		inp_distances[i] = atof(data_mx[i][1].c_str());
	}
	wutils::sort2(inp_redshifts, inp_distances);
	y2_splines.resize(nrows);
	wutils::spline(inp_redshifts, inp_distances, 1.e99, 1.e99, y2_splines);
}
