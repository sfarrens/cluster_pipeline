#include<iostream>
#include<iomanip>
#include <string.h>
#include <stdio.h>
#include "../lib/fits.hpp"          //*FITSIO packages*//

int main(int argc, char *argv[]){
  fitsfile *fptr; /*FITS file pointer, defined in fitsio.h*/
  char *val,value[1000],nullstr[]="*";
  char keyword[FLEN_KEYWORD],colname[FLEN_VALUE];
  int status = 0; /*CFITSIO status value MUST be initialized to zero!*/
  int hdunum,hdutype,ncols,anynul,dispwidth[1000],firstcol,lastcol=0,linewidth,typecode;
  long nrows,repeat,width;
  if(argc==2 || argc==3){
    if(!fits_open_file(&fptr,argv[1],READONLY,&status)){
      if(fits_get_hdu_num(fptr,&hdunum)==1) fits_movabs_hdu(fptr,2,&hdutype,&status);
      else fits_get_hdu_type(fptr,&hdutype,&status); /*Get the HDU type*/
      if (hdutype==IMAGE_HDU) std::cout<<"Error: this program only displays tables, not images\n"<<std::flush;
      else{
	fits_get_num_rows(fptr,&nrows,&status);
	fits_get_num_cols(fptr,&ncols,&status);
	while(lastcol<ncols){
	  linewidth=0;
	  firstcol=lastcol+1;
	  for(lastcol=firstcol;lastcol<=ncols;lastcol++){
	    fits_get_col_display_width(fptr,lastcol,&dispwidth[lastcol],&status);
	    linewidth+=dispwidth[lastcol]+1;
	  }
	  if(lastcol>firstcol)lastcol--; /*the last col didn't fit*/
	  std::cout<<"\t#"<<std::flush; /*print column names as column headers*/
	  for (int ii=firstcol;ii<=lastcol;ii++) {
	    fits_make_keyn("TTYPE",ii,keyword,&status);
	    fits_read_key(fptr,TSTRING,keyword,colname,NULL,&status);
	    colname[dispwidth[ii]]='\0'; /*truncate long names*/
	    std::cout<<colname<<"\t\t"<<std::flush;
	  }
	  std::cout<<"\n"<<std::flush;
	  val=value; 
	  for(int jj=1;jj<=nrows && !status;jj++){
	    for(int ii=firstcol;ii<=lastcol;ii++){
	      fits_get_coltype(fptr,ii,&typecode,&repeat,&width,&status);
	      if(typecode==16){
		fits_read_col(fptr,TSTRING,ii,jj,1,1,nullstr,&val,&anynul,&status);
		std::cout<<value<<"\t"<<std::flush;
	      }
	      else{
		if(argc==2){
		  for (int kk=1;kk<=repeat;kk++){
		    if(fits_read_col_str(fptr,ii,jj,kk,1,nullstr,&val,&anynul,&status)) break;
		    std::cout<<value<<"\t"<<std::flush;
		  }
		}
		else if(argc==3){
		  int kk=atoi(argv[2]);
		  if(kk<=repeat){
		    if(fits_read_col_str(fptr,ii,jj,kk,1,nullstr,&val,&anynul,&status)) break;
		    std::cout<<value<<"\t"<<std::flush;
		  }
		  else{
		    std::cout<<"ERROR: Vector only contains "<<repeat<<" elements\n"<<std::flush;
		    exit(-1);
		  }
		}
		else{
		  std::cout<<"ERROR: Incorrect number of arguments\n"<<std::flush;
		  exit(-1);
		}
	      }
	    }
	    std::cout<<"\n"<<std::flush;
	  }
	}
      }
      fits_close_file(fptr,&status);
    } 
    if(status)fits_report_error(stderr,status); /*print any error message*/
    return(status);
  }
  else{
    std::cout<<"Usage:  tablist filename[ext][col filter][row filter] \n"<<std::flush;
    std::cout<<"\n"<<std::flush;
    std::cout<<"List the contents of a FITS table \n"<<std::flush;
    std::cout<<"\n"<<std::flush;
    std::cout<<"Examples: \n"<<std::flush;
    std::cout<<"  tablist tab.fits[GTI]           - list the GTI extension\n"<<std::flush;
    std::cout<<"  tablist tab.fits[1][#row < 101] - list first 100 rows\n"<<std::flush;
    std::cout<<"  tablist tab.fits[1][col X;Y]    - list X and Y cols only\n"<<std::flush;
    std::cout<<"  tablist tab.fits[1][col -PI]    - list all but the PI col\n"<<std::flush;
    std::cout<<"  tablist tab.fits[1][col -PI][#row < 101]  - combined case\n"<<std::flush;
    std::cout<<"  tablist tab.fits[1][col OMAG] 3 - list element 3 from 5D column vector OMAG\n"<<std::flush;
    std::cout<<"\n"<<std::flush;
    std::cout<<"Display formats can be modified with the TDISPn keywords.\n"<<std::flush;
    return(0);
  }
}
