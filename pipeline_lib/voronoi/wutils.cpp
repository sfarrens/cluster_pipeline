/*
	Code:			wutils.h
	Author:         Walter A. Santos Jr. (walter.augusto@gmail.com)
	Based on:		Some functions based on Sam Farrens's and UCLs Weak Lensing Pipeline (WLP)
					via Filipe Abdalla
	Last update:	08/08/2012
	Version:		1.2
	Description: 	Utility codes to use on other routines, like Voronoi3D.
	Usage:			Include "wutils.h" to use in other C/C++ codes.
*/

#include "wutils.h"

std::string wutils::replace_tag(const std::string &in, const std::string &oldtag, const std::string &newtag){
	std::string file_in=in;
	std::string file_out=in;
	size_t found;

	found = file_in.find(oldtag);
	
	if(found!=0){
		file_out.replace(found, oldtag.length(), newtag);
	}
	else{ 
		file_out = wutils::tag_outfile(file_in, newtag);
	}
	return file_out;
}

void wutils::sort2(std::vector<double> &arr, std::vector<double> &brr)
{
	const int M=7,NSTACK=50;
	int i,ir,j,k,jstack=-1,l=0;
	double a,b;
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


void wutils::splint(std::vector<double> &xa, std::vector<double> &ya, std::vector<double> &y2a, const double x, double &y)
{
	int k;
	double h,b,a;

	int n=xa.size();
	int klo=0;
	int khi=n-1;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) exit(-1);
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]
		+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

void wutils::spline(std::vector<double> &x, std::vector<double> &y, const double yp1, const double ypn,
	std::vector<double> &y2)
{
	int i,k;
	double p,qn,sig,un;

	int n=y2.size();
	std::vector<double> u(n-1);
	if (yp1 > 0.99e30)
		y2[0]=u[0]=0.0;
	else {
		y2[0] = -0.5;
		u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}
	for (i=1;i<n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for (k=n-2;k>=0;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
}

std::string wutils::IntToStr(int n){
  std::ostringstream result;
  result << n;
  return result.str();
}

// given an array of double values with size num_points, find min and max values
void wutils::minmax(std::vector<double>& list, const int num_points, double& min, double& max){
	int n;

	if (num_points % 2) {
		min = list[0];
		max = list[0];
		n = 1;
	} else {
		if (list[0] > list[1]) {
			max = list[0];
			min = list[1];
		} else {
			min = list[0];
			max = list[1];
		}
		n = 2;
	}
	for (int i = n; i < num_points; i += 2) {
		if (list[i] > list[i+1]) {
			if (!(max > list[i]))
				max = list[i];
			if (!(min < list[i+1]))
				min = list[i+1];
		} else {
			if (!(max > list[i+1]))
				max = list[i+1];
			if (!(min < list[i]))
				min = list[i];
		}
	}
}



// given a multidimensional matrix of double values, find min and max boundaries
void wutils::find_boundaries(std::vector< std::vector< double > > &input_mx, const int num_dims, const int num_points, 
								std::vector<double>& min, std::vector<double>& max){

	std::vector<double> column;
	column.resize(num_points);

	for (int k=0; k < num_dims; k++){
		//copying each dimension column
		for (int i=0; i<num_points; i++)
			column[i] = input_mx[k][i];
		wutils::minmax(column, num_points, min[k], max[k]);
	}
}

// given an input filename string, return an output filename with a suffix to be appended
std::string wutils::tag_outfile(const std::string &in, const std::string &tag){
	std::string file_in=in;
	std::string file_out;
	size_t found,check;
	check = file_in.find(".");
	
	if(check!=0){
		found = file_in.find_last_of(".");
		file_out = file_in.substr(0,found);
		file_out += "." + tag;
		file_out += file_in.substr(found);
	}
	else{ 
		file_out = file_in+"."+tag;
	}
	return file_out;
}

std::string wutils::check_file_type(const std::string &fname){
  std::string c_line,r_line;
  std::size_t find_fits,find_ascii;
  c_line="file "+fname+" > file_type_check";
  system(c_line.c_str());
  std::ifstream test_file("file_type_check");
  getline(test_file,r_line);
  test_file.close();
  system("rm file_type_check");
  find_fits=r_line.find("FITS");
  find_ascii=r_line.find("ASCII");
  if(find_fits!=std::string::npos) return "fits";
  else if(find_ascii!=std::string::npos) return "ascii";
  else return "invalid";
}

void wutils::skip_comments(std::ifstream& fstr){
  char a;
  fstr>>a;
  if(a=='#'){
    fstr.ignore(2147483647,'\n');
    skip_comments(fstr);
  } else fstr.unget();
}


int wutils::count_cols(const std::string &fname){
  int n_cols;
  std::string line,temp;
  std::vector<std::string> col_count;
  
  std::ifstream test_file;
	test_file.open(fname.c_str());
  
  if (!test_file.is_open()){
		std::cout << "ERROR: Unable to open input file " << fname << " !\n";
		exit(-1);
	}
  //skip_comments(test_file);
  getline(test_file,line);
  std::stringstream line_stream(line,std::stringstream::in);
  test_file.close();
  while(!line_stream.eof()){
    line_stream>>temp;
    col_count.push_back(temp);
    line_stream>>std::ws;
  }
  n_cols=col_count.size();
  line.clear();
  temp.clear();
  line_stream.str("");
  col_count.clear();
  return n_cols;
}

void wutils::read_ascii(const std::string &fname, std::vector< std::vector< std::string > > &data_mx, int& nrows, int& ncols){
  ncols = count_cols(fname);
  
  for(int i=0;i<ncols;i++) 
	data_mx.push_back(std::vector<std::string>());
  std::string d_var;
  std::ifstream read_file(fname.c_str());
  skip_comments(read_file);
  for(;;){
    read_file>>d_var;
    if(read_file.eof()) break;
    data_mx[0].push_back(d_var);
    for(int i=1;i<ncols;i++){
      read_file>>d_var;
      data_mx[i].push_back(d_var);
    }
  }
  read_file.close();
  nrows = data_mx[0].size();
  std::cout<<"Finished reading ASCII file with "<<ncols<<" columns and "<<nrows<<" rows.\n"<<std::flush;
}


void wutils::read_fits(const std::string &fname, std::vector< std::vector< std::string > > &data_mx, int& nrows, int& ncols){
  fitsfile *fptr;
  int status=0,hdunum,hdutype,anynul;
  char *val,nullstr[]="*";
  val = new char[1000];
  long n_rows;
  if(!fits_open_file(&fptr,fname.c_str(),READONLY,&status)){
    if(fits_get_hdu_num(fptr,&hdunum)==1) fits_movabs_hdu(fptr,2,&hdutype,&status);
    else fits_get_hdu_type(fptr,&hdutype,&status);
    if(hdutype==IMAGE_HDU) std::cout<<"Error: this program only displays tables, not images.\n"<<std::flush;
    else{
      fits_get_num_rows(fptr,&n_rows,&status);
	  nrows = (int) n_rows;
      fits_get_num_cols(fptr,&ncols,&status);
      for(int i=0;i<ncols;i++) 
		data_mx.push_back(std::vector<std::string>());
        for(int i=1;i<=nrows && !status;i++){
	for(int j=1;j<=ncols;j++){
	  if(fits_read_col_str(fptr,j,i,1,1,nullstr,&val,&anynul,&status)) break;
	  data_mx[j-1].push_back(val);
	}
      }
    }
  }
  delete [] val;
  std::cout<<"Finished reading FITS file with "<<ncols<<" columns and "<<data_mx[0].size()<<" rows.\n"<<std::flush;
}



// read a file, saving the data into a string matrix nrows x ncol, Also provides nrows, ncols as int
// Base on Sam's code
void wutils::readFile(const std::string& filename, std::vector< std::vector< std::string > > &data_mx, int& nrows, int& ncols){
	std::string file_type;
	file_type = wutils::check_file_type(filename);
	
	if(!strcmp(file_type.c_str(),"ascii")) 
	read_ascii(filename, data_mx, nrows, ncols);
	else if(!strcmp(file_type.c_str(),"fits")) 
	read_fits(filename, data_mx, nrows, ncols);	
}

						
// functions to alloc/free vectors/matrices

double **wutils::alloc_2dmatrix(int rows, int cols){
	double **mx;
	mx = new double*[rows];
	for (int i=0; i < rows; i++)
		mx[i] = new double[cols];
	return mx;
}

void wutils::free_2dmatrix(double **mx, int rows){
	for (int i=0; i < rows; i++)
		delete [] mx[i];
	delete [] mx;
}

double ***wutils::alloc_3dmatrix(int x, int y, int z){
	double ***mx;

	mx = new double**[x];
	for (int i=0; i<x; ++i){
		mx[i] = new double*[y];
		for(int j=0; j<y; ++j){
			mx[i][j] = new double[z];
		}
	}
	return mx;
}

void wutils::free_3dmatrix(double ***mx, int x, int y){
	for (int i=0; i<x; ++i){
		for (int j=0; j<y; ++j)
			delete [] mx[i][j];
		delete [] mx[i];
	}
	delete [] mx;
}
