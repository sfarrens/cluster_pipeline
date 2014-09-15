//*************************************************//
//***merge_fits written by Samuel Farrens (2012)***//
//*************************************************//
//***********LAST UPDATE: 02-08-2012***************//
//*************************************************//
#include "../lib/include.hpp"       //*Include C++ packages and cosmology routines*//
#include "../lib/fits.hpp"          //*FITSIO packages*//

class do_stuff{
private:
public:  
  std::string out_file_name;
  std::vector<int> new_col;
  std::vector<std::string> in_file_name;
  void write_fits_file();
  void merge_fits(const std::string &fname,int &new_cols);
};

void do_stuff::write_fits_file(){
  fitsfile *fptr,*ofptr;
  int hdunum,hdutype,ncols,status=0,create=TRUE;
  fits_create_file(&ofptr,out_file_name.c_str(),&status); 
  fits_open_file(&fptr,in_file_name[0].c_str(),READWRITE,&status);
  if(fits_get_hdu_num(fptr,&hdunum)==1) fits_movabs_hdu(fptr,2,&hdutype,&status);
  else fits_get_hdu_type(fptr,&hdutype,&status);
  if(hdutype==IMAGE_HDU){
    std::cout<<"Error: this code only reads binary tables, not images.\n"<<std::flush;
    exit(-1);
  }
  else fits_copy_file(fptr,ofptr,1,1,1,&status);
  fits_close_file(fptr,&status);
  fits_close_file(ofptr,&status); 
}

void do_stuff::merge_fits(const std::string &fname,int &new_cols){
  fitsfile *fptr,*ofptr;
  int hdunum1,hdunum2,hdutype1,hdutype2,ncols1,ncols2,status=0,create=TRUE;
  long nrows1,nrows2;
  fits_open_file(&fptr,out_file_name.c_str(),READWRITE,&status);
  fits_open_file(&ofptr,fname.c_str(),READONLY,&status);
  if(fits_get_hdu_num(fptr,&hdunum1)==1 && fits_get_hdu_num(ofptr,&hdunum2)==1){
    fits_movabs_hdu(fptr,2,&hdutype1,&status);
    fits_movabs_hdu(ofptr,2,&hdutype2,&status);
  }
  else{
    fits_get_hdu_type(fptr,&hdutype1,&status);
    fits_get_hdu_type(ofptr,&hdutype2,&status);
  }
  if(hdutype1==IMAGE_HDU || hdutype2==IMAGE_HDU){
    std::cout<<"Error: this code only reads binary tables, not images.\n"<<std::flush;
    exit(-1);
  }
  else{  
    fits_get_num_cols(fptr,&ncols1,&status);
    fits_get_num_cols(ofptr,&ncols2,&status);
    fits_get_num_rows(fptr,&nrows1,&status);
    fits_get_num_rows(ofptr,&nrows2,&status);
    if(nrows1!=nrows2){
      std::cout<<"ERROR!: Number of rows in file 1 does not match number of rows in file 2."<<std::flush;
      exit(-1);
    }
    if(new_cols==0) new_cols=ncols1+1;
    else if(new_cols>ncols1+1){
      std::cout<<"ERROR! File only contains "<<ncols1<<" columns.\n"<<std::flush;
      exit(-1);
    }
    for(int i=1;i<=ncols2;i++){
      fits_copy_col(ofptr,fptr,i,new_cols,create,&status);
      new_cols++;
    }
  }
  fits_close_file(fptr,&status);
  fits_close_file(ofptr,&status);
}

void help();
void error1();
void error2();
void error3();
void error4();

int main (int argc, char *argv[]){
  class do_stuff do_it;
  if(argc==1) help();
  for(int i=1;i<argc;i++){
    if(argv[i][0]=='-' && std::strlen(argv[i])>1){
      switch(argv[i][1]){
      case 'h' : 
      case 'H' :
	if(std::strlen(argv[i])==2 || (std::strlen(argv[i])>2 && std::string(argv[i])=="-help")) help();
	break;
      case 'i' :
	if(std::strlen(argv[i])==2 || (std::strlen(argv[i])>2 && std::string(argv[i])=="-inp")){
	  for(int j=i+1;j<argc;j++){
	    if(argv[j][0]!='-') do_it.in_file_name.push_back(argv[j]);
	    else break;
	  }
	}
	break;
      case 'o' :
	if(std::strlen(argv[i])==2 || (std::strlen(argv[i])>2 && std::string(argv[i])=="-out")){
	  if(i!=(argc-1) && argv[i+1][0]!='-') do_it.out_file_name=argv[i+1];
	  else error3();
	  if((i+1)!=(argc-1) && argv[i+2][0]!='-') error4();
	}
	break;
      case 'c' :
	if(std::strlen(argv[i])==2 || (std::strlen(argv[i])>2 && std::string(argv[i])=="-column")){
	  for(int j=i+1;j<argc;j++){
	    if(argv[j][0]!='-') do_it.new_col.push_back(atoi(argv[j]));
	    else break;
	  }
	  if(do_it.new_col.size()!=(do_it.in_file_name.size())-1) error2();
	}
	break;
      }		
    }
  }
  if(do_it.in_file_name.size()<2) error1();
  if(do_it.out_file_name.empty()) error3();
  if(do_it.new_col.size()<1) for(int i=0;i<do_it.in_file_name.size()-1;i++) do_it.new_col.push_back(0);
  do_it.write_fits_file();
  for(int i=1;i<do_it.in_file_name.size();i++) do_it.merge_fits(do_it.in_file_name[i],do_it.new_col[i-1]);
  return 0;
}

void help(){
  std::cout<<"|-------------------------------------------------------------------------------------|\n"<<std::flush;
  std::cout<<"| This code merges the binary table data of multiple FITS files into a single output. |\n"<<std::flush;
  std::cout<<"|-------------------------------------------------------------------------------------|\n"<<std::flush;  
  std::cout<<"| OPTIONS:\t\t\t\t\t\t\t\t\t      |\n"<<std::flush;
  std::cout<<"|\t\t\t\t\t\t\t\t\t\t      |\n"<<std::flush;
  std::cout<<"|\t-h [help]\tDisplays this help page.\t\t\t\t      |\n"<<std::flush;
  std::cout<<"|\t\t\t\t\t\t\t\t\t\t      |\n"<<std::flush;
  std::cout<<"|\t-i [input]\tInput file names should be provided as arguments following    |\n"<<std::flush;
  std::cout<<"|\t\t        this option. At least two files names must be provided for    |\n"<<std::flush;
  std::cout<<"|\t\t        merging.                                                      |\n"<<std::flush;
  std::cout<<"|\t\t\t\t\t\t\t\t\t\t      |\n"<<std::flush;
  std::cout<<"|\t-o [output]\tAn output file name should be provided as an argument         |\n"<<std::flush;
  std::cout<<"|\t\t        following this option. This option must be used in            |\n"<<std::flush;
  std::cout<<"|\t\t        conjunction with -i.                                          |\n"<<std::flush;
  std::cout<<"|\t\t\t\t\t\t\t\t\t\t      |\n"<<std::flush;
  std::cout<<"|\t-c [column]\tColumn number positions should be provided as arguments       |\n"<<std::flush;
  std::cout<<"|\t\t        following this option. If this option is not used all columns |\n"<<std::flush;
  std::cout<<"|\t\t        from the input files wil be added together in the order in    |\n"<<std::flush;
  std::cout<<"|\t\t        which they are specified following the -i option. This first  |\n"<<std::flush;
  std::cout<<"|\t\t        column number will specify the position in the output file    |\n"<<std::flush;
  std::cout<<"|\t\t        where the columns from the second file will be placed. The    |\n"<<std::flush;
  std::cout<<"|\t\t        number of column positions specified should always be one     |\n"<<std::flush;
  std::cout<<"|\t\t        less than the number of input files.                          |\n"<<std::flush;
  std::cout<<"|-------------------------------------------------------------------------------------|\n"<<std::flush;
  std::cout<<"| EXAMPLES:\t\t\t\t\t\t\t\t\t      |\n"<<std::flush;
  std::cout<<"|\t\t\t\t\t\t\t\t\t\t      |\n"<<std::flush;
  std::cout<<"|\t> merge_fits -i file1.fit file2.fit -o file3.fit                              |\n"<<std::flush;
  std::cout<<"|\t\t\t\t\t\t\t\t\t\t      |\n"<<std::flush;
  std::cout<<"| This will merge all of the columns from file1.fit and file2.fit into file3.fit.     |\n"<<std::flush;
  std::cout<<"|\t\t\t\t\t\t\t\t\t\t      |\n"<<std::flush;
  std::cout<<"|\t> merge_fits -i file1.fit file2.fit -o file3.fit -c 2                         |\n"<<std::flush;
  std::cout<<"|\t\t\t\t\t\t\t\t\t\t      |\n"<<std::flush;
  std::cout<<"| This will merge all of the columns from file1.fit and file2.fit into file3.fit with |\n"<<std::flush;
  std::cout<<"| the columns from file2.fit starting in the second column position of file3.fit.     |\n"<<std::flush;
  std::cout<<"|-------------------------------------------------------------------------------------|\n"<<std::flush;  
  exit(0);
}

void error1(){
  std::cout<<"ERROR!: Must specify at least two input files for merging.\n"<<std::flush;
  std::cout<<"Use -h option to display help.\n"<<std::flush;
  exit(-1);
}

void error2(){
  std::cout<<"ERROR!: Number of column positions must be one less than number of input files.\n"<<std::flush;
  std::cout<<"Use -h option to display help.\n"<<std::flush;
  exit(-1);
}

void error3(){
  std::cout<<"ERROR!: Output file name not specified.\n"<<std::flush;
  std::cout<<"Use -h option to display help.\n"<<std::flush;
  exit(-1);
}

void error4(){
  std::cout<<"ERROR!: Only one file name can be specified for output.\n"<<std::flush;
  std::cout<<"Use -h option to display help.\n"<<std::flush;
  exit(-1);
}
