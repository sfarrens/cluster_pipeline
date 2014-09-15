/***********************************************/
/* PHOTO-Z/CLUSTER DETECTION/ANALYSIS PIPELINE */
/***********************************************/
#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<cfloat>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<limits.h>
#include<time.h>
#include<streambuf>
#include<algorithm>
#include<vector>
#include<unistd.h>
#include<sys/stat.h>
#include <unistd.h>

/***********************************************/
/*  CLASS MEMBERS                              */
/***********************************************/

class do_stuff{
private:
public:
  int run_number,block_size,num_nets,n_boxes,num_pixels,n_max_split,n_cols,flag;
  double inter_perc;
  std::string user,run,work_dir,cat_dir,temp_dir,annz_path,lephare_path,annz_train_path,annz_valid_path,
    annz_train_file,annz_valid_file,annz_net,lephare_param,lephare_context,lephare_pdz,annz_itt,pipe_path,id_col,cf_file_list,
    cf_file_path,op_use_z,op_dir_set_up,op_make_initial_cats,op_make_lephare,op_make_annz_train,op_make_annz_test,op_make_full_cats,
    op_run_xyz,op_run_find_boundaries,op_run_main_split,op_run_split,op_run_voronoi,op_run_cleaup,log_file_name,op_spec_z,rm_temp,
    op_run_pipeline_linking,op_density_cuts,density_cut_value,op_run_find_boundaries_radec,op_run_split_radec;
  std::vector<int> jobid_block,neg_block,net_block,tempjobs;
  std::vector<std::string> cf_file_name,wts,rand_list,no_mag_string,new_mag_string,mag_filt_lim,pos_cols,mag_cols,mag_err_cols;
  void read_params(char *file_named);
  void read_file_list(const std::string &fname1,const std::string &fname2);
  void split(const std::string &str,std::vector<std::string> &tokens,const std::string &delimiter);
  void whileEOF(std::ifstream &inputfile,std::vector<std::string> &header,std::vector<std::string> &values,
		const std::string &comment_str);
  void skip_comments(std::ifstream& fstr);
  void count_cols(const std::string &fname);
  void set_up(int x);
  void dir_set_up();
  void set_up_log();
  void read_file_list();
  void pbs_initial_cats(std::string name,int i);
  void make_initial_cats(int i,int j);
  void pbs_annz_train(std::string name,int j);
  void make_annz_train(int i,int j);
  void pbs_annz_test(std::string name,int i);
  void make_annz_test(int i,int j);
  void pbs_lephare(std::string name,int i);
  void make_lephare(int i,int j);
  void pbs_full_cats(std::string name,int i);
  void make_full_cats(int i,int j);
  void pbs_xyz(std::string name,int i);
  void run_xyz(int i,int j);
  void pbs_find_boundaries(std::string name);
  void run_find_boundaries();
  void pbs_main_split(std::string name);
  void run_main_split();
  void pbs_split(std::string name,int i,int o);
  void run_split(int i,int j,int o);
  void pbs_voronoi(std::string name,int i);
  void run_voronoi(int i,int j);
  void pbs_cleanup(std::string name,int i,int k);
  void run_cleanup(int i,int j,int k);
  void pbs_pipeline_linking(std::string name,int i);
  void run_pipeline_linking(int i,int j);
  void pbs_density_cuts(std::string name,int i);
  void run_density_cuts(int i,int j);
  void temp_jobs();
  void pbs_split_radec(std::string name,int i,int o);
  void run_split_radec(int i,int j,int o);
};

/***********************************************/
/*  SPLIT FILE INTO COLUMNS                    */
/***********************************************/

void do_stuff::split(const std::string &str,std::vector<std::string> &tokens,const std::string &delimiter){
  std::string::size_type lastPos,pos;      
  lastPos=str.find_first_not_of(delimiter,0); /* pos = find first "non-delimiter" */
  pos=str.find_first_of(delimiter,lastPos);
  while (std::string::npos!=pos || std::string::npos!=lastPos){
    tokens.push_back(str.substr(lastPos, pos-lastPos)); /* Found a token, add it to the vector. */     
    lastPos=str.find_first_not_of(delimiter, pos); /* Skip delimiters. Note the "not_of" */         
    pos=str.find_first_of(delimiter, lastPos); /* Find next "non-delimiter" */            
  }
}
  
void do_stuff::whileEOF(std::ifstream &inputfile,std::vector<std::string> &header,std::vector<std::string> &values,
			const std::string &comment_str){
  std::string line;
  int hash;
  while(!inputfile.eof()){ /* While not the end of the file */
    std::getline(inputfile,line); /* Read each line */
    if(line.length() >= 1){ /* For lines with length > 1 (i.e., skip lines with no length = empty lines) */
      int hash=line.find(comment_str); /* Find lines that start with a # - indicating header information */
      if(hash!=std::string::npos) /* If you find a line starting with comment_str */
	header.push_back(line); /* Read the headers here and push to 'header' container */
      else
	values.push_back(line); /* Read the values here and push to 'values' container */
    }
  }
}

/***********************************************/
/*  SKIP LINES STARTING WITH #                 */
/***********************************************/

void do_stuff::skip_comments(std::ifstream& fstr){
  char a;
  fstr>>a;
  if(a=='#'){
    fstr.ignore(2147483647,'\n');
    skip_comments(fstr);
  } else fstr.unget();
}

/***********************************************/
/*  COUNT FILE COLUMNS                         */
/***********************************************/

void do_stuff::count_cols(const std::string &fname){
  std::string line,temp;
  std::vector<std::string> col_count;
  std::ifstream test_file(fname.c_str());
  skip_comments(test_file);
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
}

/***********************************************/
/*  READ PARAMETER FILE                        */
/***********************************************/

void do_stuff::read_params(char *file_named){
  std::ifstream file(file_named);
  std::string arch_line,up_lim_line,low_lim_line,start_line;
  std::vector<std::string> header,lines,values;
  whileEOF(file,header,lines,"#");
  file.close();
  for(int i=0;i<lines.size();++i){
    values.clear();
    split(lines[i],values," "); 
    if(values[0]=="user") user=values[1];
    if(values[0]=="run_number") run_number=atoi(values[1].c_str());
    if(values[0]=="work_dir") work_dir=values[1];
    if(values[0]=="cat_dir") cat_dir=values[1];
    if(values[0]=="block_size") block_size=atoi(values[1].c_str());
    if(values[0]=="num_nets") num_nets=atoi(values[1].c_str());
    if(values[0]=="no_mag_string") for(int j=1;j<values.size();j++) no_mag_string.push_back(values[1]);
    if(values[0]=="new_mag_string") for(int j=1;j<values.size();j++) new_mag_string.push_back(values[1]);
    if(values[0]=="mag_filt_lim") for(int j=1;j<values.size();j++) mag_filt_lim.push_back(values[j]);
    if(values[0]=="annz_itt") annz_itt=values[1];
    if(values[0]=="annz_train_path") annz_train_path=values[1];
    if(values[0]=="annz_valid_path") annz_valid_path=values[1];
    if(values[0]=="annz_train_file") annz_train_file=values[1];
    if(values[0]=="annz_valid_file") annz_valid_file=values[1];
    if(values[0]=="annz_net") annz_net=values[1];
    if(values[0]=="annz_path") annz_path=values[1];
    if(values[0]=="lephare_path") lephare_path=values[1];
    if(values[0]=="lephare_param") lephare_param=values[1];
    if(values[0]=="lephare_context") lephare_context=values[1];
    if(values[0]=="lephare_pdz") lephare_pdz=values[1];
    if(values[0]=="pipe_path") pipe_path=values[1];
    if(values[0]=="n_max_split") n_max_split=atoi(values[1].c_str());
    if(values[0]=="cf_file_path") cf_file_path=values[1];
    if(values[0]=="cf_file_list") cf_file_list=values[1];
    if(values[0]=="pos_cols") for(int j=1;j<values.size();j++) pos_cols.push_back(values[j]);
    if(values[0]=="mag_cols") for(int j=1;j<values.size();j++) mag_cols.push_back(values[j]);
    if(values[0]=="mag_err_cols") for(int j=1;j<values.size();j++) mag_err_cols.push_back(values[j]);
    if(values[0]=="density_cut_value") density_cut_value=values[1];
    if(values[0]=="num_pixels") num_pixels=atoi(values[1].c_str());
    if(values[0]=="inter_perc") inter_perc=atof(values[1].c_str());
    if(values[0]=="op_use_z") op_use_z=values[1];
    if(values[0]=="op_spec_z") op_spec_z=values[1];
    if(values[0]=="rm_temp") rm_temp=values[1];
    if(values[0]=="op_dir_set_up") op_dir_set_up=values[1];
    if(values[0]=="op_make_initial_cats") op_make_initial_cats=values[1];
    if(values[0]=="op_make_lephare") op_make_lephare=values[1];
    if(values[0]=="op_make_annz_train") op_make_annz_train=values[1];
    if(values[0]=="op_make_annz_test") op_make_annz_test=values[1];
    if(values[0]=="op_make_full_cats") op_make_full_cats=values[1];
    if(values[0]=="op_run_xyz") op_run_xyz=values[1];
    if(values[0]=="op_run_find_boundaries") op_run_find_boundaries=values[1];
    if(values[0]=="op_run_main_split") op_run_main_split=values[1];
    if(values[0]=="op_run_split") op_run_split=values[1];
    if(values[0]=="op_run_voronoi") op_run_voronoi=values[1];
    if(values[0]=="op_run_cleaup") op_run_cleaup=values[1];
    if(values[0]=="op_run_pipeline_linking") op_run_pipeline_linking=values[1];
    if(values[0]=="op_density_cuts") op_density_cuts=values[1];
    if(values[0]=="op_run_find_boundaries_radec") op_run_find_boundaries_radec=values[1];
    if(values[0]=="op_run_split_radec") op_run_split_radec=values[1];
  }
  if(strcmp(op_use_z.c_str(),"annz") && strcmp(op_use_z.c_str(),"lephare") && strcmp(op_use_z.c_str(),"incz")){
    std::cout<<op_use_z<<" IS NOT A VALID OPTION FOR <op_use_z>. OPTIONS ARE: annz lephare incz\n"<<std::flush;
    std::cout<<"PLEASE EDIT PARAMETER FILE.\n";
    exit(-1);
  }
}

/***********************************************/
/*  GENERAL SET UP FEATURES                    */
/***********************************************/

void do_stuff::set_up(int x){
  std::ostringstream run_stream,rand_stream;
  run_stream<<std::setw(4)<<std::setfill('0')<<x;
  run=run_stream.str(); 
  temp_dir=user+"_PIPELINE_RUN_"+run+"_TEMP/";
  for(int i=0;i<block_size;i++){
    jobid_block.push_back(-1);
    neg_block.push_back(-1);
  }
  for(int i=0;i<num_nets;i++)
    net_block.push_back(-1);
  for(int i=0;i<100000;i++){
    rand_stream<<i+1;
    rand_list.push_back(rand_stream.str());
    rand_stream.str("");
  }
  srand(time(NULL));
  random_shuffle(rand_list.begin(),rand_list.end());
  flag=0;
}

/***********************************************/
/*  SET UP RUN DIRECTORIES                     */
/***********************************************/

void do_stuff::dir_set_up(){
  std::string temp_str;
  temp_str=cat_dir+"PIPELINE_RUN_"+run;
  mkdir(temp_str.c_str(),S_IRWXU|S_IRGRP|S_IXGRP);
  temp_str=cat_dir+"PIPELINE_RUN_"+run+"/NEW_CATALOGUES";
  mkdir(temp_str.c_str(),S_IRWXU|S_IRGRP|S_IXGRP);
  temp_str=cat_dir+"PIPELINE_RUN_"+run+"/PHOTOZ";
  mkdir(temp_str.c_str(),S_IRWXU|S_IRGRP|S_IXGRP);
  temp_str=cat_dir+"PIPELINE_RUN_"+run+"/PHOTOZ/NETWORKS";
  mkdir(temp_str.c_str(),S_IRWXU|S_IRGRP|S_IXGRP);
  temp_str=cat_dir+"PIPELINE_RUN_"+run+"/PHOTOZ/PDZ";
  mkdir(temp_str.c_str(),S_IRWXU|S_IRGRP|S_IXGRP);
  temp_str=cat_dir+"PIPELINE_RUN_"+run+"/MAGS";
  mkdir(temp_str.c_str(),S_IRWXU|S_IRGRP|S_IXGRP);
  temp_str=cat_dir+"PIPELINE_RUN_"+run+"/VORONOI";
  mkdir(temp_str.c_str(),S_IRWXU|S_IRGRP|S_IXGRP);
  temp_str=cat_dir+"PIPELINE_RUN_"+run+"/VORONOI/SPLITS";
  mkdir(temp_str.c_str(),S_IRWXU|S_IRGRP|S_IXGRP);
  temp_str=cat_dir+"PIPELINE_RUN_"+run+"/VORONOI/SPLITS_RADEC";
  mkdir(temp_str.c_str(),S_IRWXU|S_IRGRP|S_IXGRP);
  temp_str=cat_dir+"PIPELINE_RUN_"+run+"/FOF";
  mkdir(temp_str.c_str(),S_IRWXU|S_IRGRP|S_IXGRP);
  temp_str=cat_dir+"PIPELINE_RUN_"+run+"/SCRIPTS";
  mkdir(temp_str.c_str(),S_IRWXU|S_IRGRP|S_IXGRP);
  temp_str=cat_dir+"PIPELINE_RUN_"+run+"/SCRIPTS/SPLITS";
  mkdir(temp_str.c_str(),S_IRWXU|S_IRGRP|S_IXGRP);
  temp_str=cat_dir+"PIPELINE_RUN_"+run+"/SCRIPTS/SPLITS_RADEC";
  mkdir(temp_str.c_str(),S_IRWXU|S_IRGRP|S_IXGRP);
  temp_str=cat_dir+"PIPELINE_RUN_"+run+"/LOG";
  mkdir(temp_str.c_str(),S_IRWXU|S_IRGRP|S_IXGRP);
}

/***********************************************/
/*  SET UP PIPELINE LOG                        */
/***********************************************/

void do_stuff::set_up_log(){
  struct tm *current;
  time_t now;
  time(&now);
  current=localtime(&now);
  log_file_name=cat_dir+"PIPELINE_RUN_"+run+"/LOG/log_file.txt";
  std::ofstream outfile;
  outfile.open(log_file_name.c_str());
  outfile<<"#----------------------------------------------------------------------------------------------#\n"<<std::flush;
  outfile<<"#--- PIPELINE LOG FILE CREATED "<<std::setw(2)<<std::setfill('0')<<current->tm_mday<<"/"<<std::setw(2)
	 <<std::setfill('0')<<current->tm_mon+1<<"/"<<std::setw(2)<<std::setfill('0')<<current->tm_year+1900
	 <<" AT "<<std::setw(2)<<std::setfill('0')<<current->tm_hour<<":"<<std::setw(2)<<std::setfill('0')<<current->tm_min
	 <<":"<<std::setw(2)<<std::setfill('0')<<current->tm_sec<<" -----------------------------------------#\n"<<std::flush;
  outfile<<"#----------------------------------------------------------------------------------------------#\n\n"<<std::flush;
  outfile<<"- USER NAME: "<<user<<"\n"<<std::flush;
  outfile<<"- RUN NUMBER: "<<run<<"\n"<<std::flush; 
  outfile.close();
}

/***********************************************/
/*  READ LIST OF CATALOGUE FILES               */
/***********************************************/

void do_stuff::read_file_list(){
  std::string s_var;
  std::ifstream read_file(cf_file_list.c_str());
  if(read_file.fail()){
    std::ofstream outlogfile;
    outlogfile.open(log_file_name.c_str(),std::ios_base::app);
    outlogfile<<"- ERROR! [Failed to open file: "<<cf_file_list<<"]\n"<<std::flush;
    outlogfile<<"- PIPELINE ABORTED...\n"<<std::flush;
    outlogfile.close();
    std::cout<<"- ERROR! [Failed to open file: "<<cf_file_list<<"]\n"<<std::flush;
    std::cout<<"- PIPELINE ABORTED...\n"<<std::flush;
    exit(-1);
  }
  for(;;){
    read_file>>s_var;
    if(read_file.eof()) break;
    cf_file_name.push_back(s_var);
  }
  read_file.close();
  std::ofstream outlogfile;
  outlogfile.open(log_file_name.c_str(),std::ios_base::app);
  outlogfile<<"- PREPARING PIPELINE FOR ["<<cf_file_name.size()<<"] CATALOGUE PIECES\n"<<std::flush;
  outlogfile.close();
}

/***********************************************/
/*  SET UP ANNZ AND LEPHARE INPUT FILES        */
/***********************************************/

/*Generate PBS Script*/
void do_stuff::pbs_initial_cats(std::string name, int i){
  std::ofstream res(name.c_str());
  res<<"#!/bin/tcsh\n";
  res<<"#PBS -lwalltime=4:00:00\n";
  res<<"#PBS -lnodes=1:ppn=1\n";
  res<<"#echo $PBS_JOBID\n";
  res<<"cd "+work_dir+"\n";
  res<<"mkdir "+temp_dir+"\n";
  std::ostringstream mag_col_stream,mag_err_col_stream,mag_col_stream2,mag_err_col_stream2,no_mag_stream,
    new_mag_stream,mag_lim_stream;
  for(int k=0;k<mag_cols.size();k++){
    mag_col_stream<<"$"<<mag_cols[k];
    mag_err_col_stream<<"$"<<mag_err_cols[k];
    mag_col_stream2<<mag_cols[k]<<" ";
    mag_err_col_stream2<<mag_err_cols[k]<<" ";
    if(k!=mag_cols.size()-1){
      mag_col_stream<<",";
      mag_err_col_stream<<",";
    }
  }
  for(int k=0;k<no_mag_string.size();k++) no_mag_stream<<no_mag_string[k]<<" ";
  for(int k=0;k<new_mag_string.size();k++) new_mag_stream<<new_mag_string[k]<<" ";
  for(int k=0;k<mag_filt_lim.size();k++) mag_lim_stream<<mag_filt_lim[k]<<" ";
  res<<"cp "+cf_file_path+cf_file_name[i]+" "+temp_dir+"\n";
  res<<pipe_path+"mag_replace -i "+temp_dir+cf_file_name[i]+" -s "+no_mag_stream.str()+"-r "+mag_lim_stream.str()+"-c "
    +mag_col_stream2.str()+mag_err_col_stream2.str()+"> "+temp_dir+cf_file_name[i]+".mags_annz_clean\n";
  res<<pipe_path+"mag_replace -i "+temp_dir+cf_file_name[i]+" -s "+no_mag_stream.str()+"-r "+new_mag_stream.str()+"-c "
    +mag_col_stream2.str()+mag_err_col_stream2.str()+"> "+temp_dir+cf_file_name[i]+".mags_lephare_clean\n";
  if(!strcmp(op_spec_z.c_str(),"yes")){
    res<<"awk 'NR>1 {print "+mag_col_stream.str()+","+mag_err_col_stream.str()+",$"+pos_cols[3]+"}' "
      +temp_dir+cf_file_name[i]+".mags_annz_clean > "+temp_dir+cf_file_name[i]+".mags_annz\n";
    res<<"awk 'NR>1 {print $"+pos_cols[0]+","+mag_col_stream.str()+","+mag_err_col_stream.str()+","+lephare_context
      +",$"+pos_cols[3]+"}' "+temp_dir+cf_file_name[i]+".mags_lephare_clean > "+temp_dir+cf_file_name[i]+".mags_lephare\n";
  }
  else{
    res<<"awk 'NR>1 {print "+mag_col_stream.str()+","+mag_err_col_stream.str()+"}' "+temp_dir+cf_file_name[i]
      +".mags_annz_clean > "+temp_dir+cf_file_name[i]+".mags_annz\n";
    res<<"awk 'NR>1 {print $"+pos_cols[0]+","+mag_col_stream.str()+","+mag_err_col_stream.str()+"}' "
      +temp_dir+cf_file_name[i]+".mags_lephare_clean > "+temp_dir+cf_file_name[i]+".mags_lephare\n";
  }
  res<<"cp "+temp_dir+cf_file_name[i]+".mags_annz "+temp_dir+cf_file_name[i]+".mags_lephare "+cat_dir
    +"PIPELINE_RUN_"+run+"/MAGS/\n";
  res<<"rm "+temp_dir+cf_file_name[i]+".mags_annz_clean "+temp_dir+cf_file_name[i]+".mags_lephare_clean\n";
  if(!strcmp(rm_temp.c_str(),"yes"))
    res<<"rm "+temp_dir+cf_file_name[i]+" "+temp_dir+cf_file_name[i]+".mags_annz "+temp_dir+cf_file_name[i]+".mags_lephare\n";
  res.close();
}

/*Submit Jobs*/
void do_stuff::make_initial_cats(int i,int j){
  char buff[512];
  std::ostringstream id_stream,idpart;
  std::string name1,name2,command,result;
  name1=cf_file_path+cf_file_name[i];
  std::ifstream read_file(name1.c_str());
  if(read_file.fail()){
    std::ofstream outlogfile;
    outlogfile.open(log_file_name.c_str(),std::ios_base::app);
    outlogfile<<"- ERROR! [Failed to open file: "<<cf_file_path<<cf_file_name[i]<<"]\n"<<std::flush;
    outlogfile<<"- PIPELINE ABORTED...\n"<<std::flush;
    outlogfile.close();
    std::cout<<"- ERROR! [Failed to open file: "<<cf_file_path<<cf_file_name[i]<<"]\n"<<std::flush;
    std::cout<<"- PIPELINE ABORTED...\n"<<std::flush;
    exit(-1);
  }
  read_file.close();
  name1=cat_dir+"PIPELINE_RUN_"+run+"/MAGS/"+cf_file_name[i]+".mags_annz";
  name2=cat_dir+"PIPELINE_RUN_"+run+"/MAGS/"+cf_file_name[i]+".mags_lephare";
  std::ifstream test_file1(name1.c_str());
  std::ifstream test_file2(name2.c_str());
  if(test_file1.fail() && test_file2.fail()){
    name1=cat_dir+"PIPELINE_RUN_"+run+"/SCRIPTS/make_initial_cats."+cf_file_name[i]+".sh";
    pbs_initial_cats(name1,i);
    idpart<<jobid_block[j];
    if(jobid_block[j]==-1) command="qsub "+name1;
    else command="qsub -W depend=afterok:"+idpart.str()+" "+name1;
    FILE *in=popen(command.c_str(),"r");
    fgets(buff, sizeof(buff), in);
    result=buff;
    for(int k=0;k<result.size();k++){
      if(isdigit(result[k]))
	id_stream<<result[k];
    }
    jobid_block[j]=atoi(id_stream.str().c_str());
    pclose(in);
    std::ofstream outlogfile;
    outlogfile.open(log_file_name.c_str(),std::ios_base::app);
    outlogfile<<">> JOB["<<jobid_block[j]<<"]: "<<command<<"\n"<<std::flush;
    outlogfile.close();
    command.clear();
    result.clear();
    idpart.str("");
    id_stream.str("");
  }
  test_file1.close();
  test_file2.close();
  flag=1;
}

/***********************************************/
/*  TRAIN ANNZ NETWORK                         */
/***********************************************/

/*Generate PBS Script*/
void do_stuff::pbs_annz_train(std::string name,int j){
  std::string test_file;
  std::ofstream res(name.c_str());
  res<<"#!/bin/tcsh\n";
  res<<"#PBS -lwalltime=99:00:00\n";
  res<<"#PBS -lnodes=1:ppn=1\n";
  res<<"#echo $PBS_JOBID\n";
  res<<"cd "+work_dir+"\n";
  res<<"mkdir "+temp_dir+"\n";
  res<<"if (! -e "+temp_dir+annz_train_file+") cp "+annz_train_path+annz_train_file+" "+temp_dir+"\n";
  res<<"if (! -e "+temp_dir+annz_valid_file+") cp "+annz_valid_path+annz_valid_file+" "+temp_dir+"\n";
  res<<"echo "+annz_itt+" 0 | "+annz_path+"annz_train "+annz_net+" "+temp_dir+annz_train_file+" "+temp_dir
    +annz_valid_file+" "+temp_dir+annz_train_file+".wts."+wts[j]+" "+rand_list[j]+"\n";
  res<<"cp "+temp_dir+annz_train_file+".wts."+wts[j]+" "+cat_dir+"PIPELINE_RUN_"+run+"/PHOTOZ/NETWORKS/\n";
  if(!strcmp(rm_temp.c_str(),"yes"))
    res<<"rm "+temp_dir+annz_train_file+" "+temp_dir+annz_valid_file+" "+temp_dir+annz_train_file+".wts."+wts[j]+"\n";
  res.close();
}

/*Submit Jobs*/
void do_stuff::make_annz_train(int i,int j){
  char buff[512];
  std::ostringstream wts_stream,idpart,id_stream;
  std::string name1,name2,name3,result,command;
  name2=annz_train_path+annz_train_file;
  name3=annz_valid_path+annz_valid_file;
  std::ifstream read_file2(name2.c_str());
  std::ifstream read_file3(name3.c_str());
  if(read_file2.fail()){
    std::ofstream outlogfile;
    outlogfile.open(log_file_name.c_str(),std::ios_base::app);
    outlogfile<<"- ERROR! [Failed to open file: "<<name2<<"]\n"<<std::flush;
    outlogfile<<"- PIPELINE ABORTED...\n"<<std::flush;
    outlogfile.close();
    std::cout<<"- ERROR! [Failed to open file: "<<name2<<"]\n"<<std::flush;
    std::cout<<"- PIPELINE ABORTED...\n"<<std::flush;
    exit(-1);
  }
  if(read_file3.fail()){
    std::ofstream outlogfile;
    outlogfile.open(log_file_name.c_str(),std::ios_base::app);
    outlogfile<<"- ERROR! [Failed to open file: "<<name3<<"]\n"<<std::flush;
    outlogfile<<"- PIPELINE ABORTED...\n"<<std::flush;
    outlogfile.close();
    std::cout<<"- ERROR! [Failed to open file: "<<name3<<"]\n"<<std::flush;
    std::cout<<"- PIPELINE ABORTED...\n"<<std::flush;
    exit(-1);
  }
  read_file2.close();
  read_file3.close();
  wts_stream<<std::setw(4)<<std::setfill('0')<<i+1;
  wts.push_back(wts_stream.str());
  name1=cat_dir+"PIPELINE_RUN_"+run+"/PHOTOZ/NETWORKS/"+annz_train_file+".wts."+wts[i]+"";
  std::ifstream test_file(name1.c_str());
  if(test_file.fail()){ 
    name1=cat_dir+"PIPELINE_RUN_"+run+"/SCRIPTS/make_annz_train."+annz_train_file+"."+wts_stream.str()+".sh";
    pbs_annz_train(name1,i);
    idpart<<jobid_block[j];
    if(jobid_block[j]==-1) command="qsub "+name1;
    else command="qsub -W depend=afterok:"+idpart.str()+" "+name1;      
    FILE *in=popen(command.c_str(),"r");
    fgets(buff, sizeof(buff), in);
    result=buff;
    for(int l=0;l<result.size();l++){
      if(isdigit(result[l]))
	id_stream<<result[l];
    }
    jobid_block[j]=atoi(id_stream.str().c_str());
    net_block[j]=atoi(id_stream.str().c_str());
    pclose(in);
    std::ofstream outlogfile;
    outlogfile.open(log_file_name.c_str(),std::ios_base::app);
    outlogfile<<">> JOB["<<jobid_block[j]<<"]: "<<command<<"\n"<<std::flush;
    outlogfile.close();
    command.clear();
    result.clear();
    wts_stream.str("");
    idpart.str("");
    id_stream.str("");
  }
  test_file.close();
  flag=1;
}

/***********************************************/
/*  RUN ANNZ                                   */
/***********************************************/

/*Generate PBS Script*/
void do_stuff::pbs_annz_test(std::string name,int i){
  std::string test_file;
  std::ofstream res(name.c_str());
  res<<"#!/bin/tcsh\n";
  res<<"#PBS -lwalltime=99:00:00\n";
  res<<"#PBS -lnodes=1:ppn=1\n";
  res<<"#echo $PBS_JOBID\n";
  res<<"cd "+work_dir+"\n";
  res<<"mkdir "+temp_dir+"\n";
  for(int j=1;j<=num_nets;j++){ 
    res<<"if (! -e "+temp_dir+annz_train_file+".wts."+wts[j-1]+") cp "+cat_dir+"PIPELINE_RUN_"
      +run+"/PHOTOZ/NETWORKS/"+annz_train_file+".wts."+wts[j-1]+" "+temp_dir+"\n";
  }
  if(!strcmp(op_make_initial_cats.c_str(),"no"))
    res<<"cp "+cf_file_path+cf_file_name[i]+" "+temp_dir+cf_file_name[i]+".mags_annz\n";
  else
    res<<"if (! -e "+temp_dir+cf_file_name[i]+".mags_annz) cp "+cat_dir+"PIPELINE_RUN_"+run+"/MAGS/"
      +cf_file_name[i]+".mags_annz "+temp_dir+"\n";
  res<<annz_path+"annz_test "+temp_dir+cf_file_name[i]+".mags_annz "+temp_dir+cf_file_name[i]+".annz";
  for(int j=1;j<=num_nets;j++){ 
    res<<" "+temp_dir+annz_train_file+".wts."+wts[j-1];
  }
  res<<"\n";
  res<<"cp "+temp_dir+cf_file_name[i]+".annz "+cat_dir+"PIPELINE_RUN_"+run+"/PHOTOZ/\n";
  if(!strcmp(rm_temp.c_str(),"yes")){
    res<<"rm "+temp_dir+cf_file_name[i]+".mags_annz "+temp_dir+cf_file_name[i]+".annz";
    for(int j=1;j<=num_nets;j++){ 
      res<<" "+temp_dir+annz_train_file+".wts."+wts[j-1];
    }
    res<<"\n";
  }
  res.close();
}

/*Submit Jobs*/
void do_stuff::make_annz_test(int i,int j){
  char buff[512];
  int net_flag;
  std::ostringstream netpart,idpart,id_stream;
  std::string name,result,command;
  name=cat_dir+"PIPELINE_RUN_"+run+"/PHOTOZ/"+cf_file_name[i]+".annz";
  std::ifstream test_file(name.c_str());
  for(int m=0;m<num_nets;m++){
    if(net_block[m]>-1)
      netpart<<":"<<net_block[m];
  }
  if(test_file.fail()){
    name=cat_dir+"PIPELINE_RUN_"+run+"/SCRIPTS/make_annz_test."+cf_file_name[i]+".sh";
    pbs_annz_test(name,i);
    idpart<<jobid_block[j];
    net_flag=0;
    for(int pp=0;pp<num_nets;pp++){
      if(jobid_block[j]==net_block[pp])
	net_flag++;
    }
    if(jobid_block==neg_block) command="qsub "+name;
    else if(net_flag==0)
      command="qsub -W depend=afterok:"+idpart.str()+netpart.str()+" "+name;  
    else
      command="qsub -W depend=afterok"+netpart.str()+" "+name;  
    FILE *in=popen(command.c_str(),"r");
    fgets(buff, sizeof(buff), in);
    result=buff;
    for(int l=0;l<result.size();l++){
      if(isdigit(result[l]))
	id_stream<<result[l];
    }
    jobid_block[j]=atoi(id_stream.str().c_str());
    pclose(in);
    std::ofstream outlogfile;
    outlogfile.open(log_file_name.c_str(),std::ios_base::app);
    outlogfile<<">> JOB["<<jobid_block[j]<<"]: "<<command<<"\n"<<std::flush;
    outlogfile.close();
    command.clear();
    result.clear();
    idpart.str("");
    id_stream.str("");
  }
  test_file.close();
  flag=1;
}

/***********************************************/
/*  RUN LEPHARE                                */
/***********************************************/

/*Generate PBS Script*/
void do_stuff::pbs_lephare(std::string name, int i){
  std::string test_file;
  std::ofstream res(name.c_str());
  res<<"#!/bin/tcsh\n";
  res<<"#PBS -lwalltime=99:00:00\n";
  res<<"#PBS -lnodes=1:ppn=1\n";
  res<<"#echo $PBS_JOBID\n";
  res<<"cd "+work_dir+"\n";
  res<<"mkdir "+temp_dir+"\n";
  res<<"setenv LEPHAREDIR "+lephare_path+"lephare_dev\n";
  res<<"setenv LEPHAREWORK "+lephare_path+"lephare_work\n";
  res<<lephare_path+"lephare_dev/source/sedtolib -t G -c "+lephare_param+"\n";
  res<<lephare_path+"lephare_dev/source/filter -c "+lephare_param+"\n";
  res<<lephare_path+"lephare_dev/source/mag_gal -t G -c "+lephare_param+"\n";
  if(!strcmp(op_make_initial_cats.c_str(),"no"))
    res<<"cp "+cf_file_path+cf_file_name[i]+" "+temp_dir+cf_file_name[i]+".mags_lephare\n";
  else
    res<<"if (! -e "+temp_dir+cf_file_name[i]+".mags_lephare) cp "+cat_dir+"PIPELINE_RUN_"
      +run+"/MAGS/"+cf_file_name[i]+".mags_lephare "+temp_dir+"\n";
  res<<lephare_path+"lephare_dev/source/zphota -c "+lephare_param+" -CAT_IN "+temp_dir
    +cf_file_name[i]+".mags_lephare -CAT_OUT "+temp_dir+cf_file_name[i]+".lephare";
  if(!strcmp(op_spec_z.c_str(),"yes")) 
    res<<" -CAT_TYPE LONG";
  else 
    res<<" -CAT_TYPE SHORT";
  if(!strcmp(lephare_pdz.c_str(),"yes")) 
    res<<" -PDZ_OUT "+temp_dir+cf_file_name[i]+".lephare_pdz";
  res<<"\n";
  res<<"cp "+temp_dir+cf_file_name[i]+".lephare "+cat_dir+"PIPELINE_RUN_"+run+"/PHOTOZ/\n";
  if(!strcmp(lephare_pdz.c_str(),"yes")){
    res<<"cp "+temp_dir+cf_file_name[i]+".lephare_pdz.pdz "+temp_dir+cf_file_name[i]+".lephare_pdz.zph "+temp_dir+cf_file_name[i]
      +".lephare_pdz.mod "+temp_dir+cf_file_name[i]+".lephare_pdz.abs02 "+temp_dir+cf_file_name[i]+".lephare_pdz.abs10 "
      +temp_dir+cf_file_name[i]+".lephare_pdz.abs14 "+cat_dir+"PIPELINE_RUN_"+run+"/PHOTOZ/PDZ/\n";
    if(!strcmp(rm_temp.c_str(),"yes"))
      res<<"rm "+temp_dir+cf_file_name[i]+".lephare_pdz.pdz "+temp_dir+cf_file_name[i]+".lephare_pdz.zph "+temp_dir+cf_file_name[i]
	+".lephare_pdz.mod "+temp_dir+cf_file_name[i]+".lephare_pdz.abs02 "+temp_dir+cf_file_name[i]+".lephare_pdz.abs10 "
	+temp_dir+cf_file_name[i]+".lephare_pdz.abs14\n";
  }
  if(!strcmp(rm_temp.c_str(),"yes"))
    res<<"rm "+temp_dir+cf_file_name[i]+".mags_lephare "+temp_dir+cf_file_name[i]+".lephare\n";
  res.close();
}

/*Submit Jobs*/
void do_stuff::make_lephare(int i,int j){
  char buff[512];
  std::ostringstream idpart,id_stream;
  std::string name,result,command;
  std::ifstream read_file1(lephare_param.c_str());
  if(read_file1.fail() && flag==0){
    std::ofstream outlogfile;
    outlogfile.open(log_file_name.c_str(),std::ios_base::app);
    outlogfile<<"- ERROR! [Failed to open file: "<<lephare_param<<"]\n"<<std::flush;
    outlogfile<<"- PIPELINE ABORTED...\n"<<std::flush;
    outlogfile.close();
    std::cout<<"- ERROR! [Failed to open file: "<<lephare_param<<"]\n"<<std::flush;
    std::cout<<"- PIPELINE ABORTED...\n"<<std::flush;
    exit(-1);
  }
  read_file1.close();
  name=cat_dir+"PIPELINE_RUN_"+run+"/PHOTOZ/"+cf_file_name[i]+".lephare";
  std::ifstream test_file(name.c_str());
  if(test_file.fail()){
    name=cat_dir+"PIPELINE_RUN_"+run+"/SCRIPTS/make_lephare."+cf_file_name[i]+".sh";
    pbs_lephare(name,i);
    idpart<<jobid_block[j];
    if(jobid_block[j]==-1) command="qsub "+name;
    else command="qsub -W depend=afterok:"+idpart.str()+" "+name;
    FILE *in=popen(command.c_str(),"r");
    fgets(buff, sizeof(buff), in);
    result=buff;
    for(int k=0;k<result.size();k++){
      if(isdigit(result[k]))
	id_stream<<result[k];
    }
    jobid_block[j]=atoi(id_stream.str().c_str());
    pclose(in);
    std::ofstream outlogfile;
    outlogfile.open(log_file_name.c_str(),std::ios_base::app);
    outlogfile<<">> JOB["<<jobid_block[j]<<"]: "<<command<<"\n"<<std::flush;
    outlogfile.close();
    command.clear();
    result.clear();
    idpart.str("");
    id_stream.str("");
  }
  test_file.close();
  flag=1;
}

/***********************************************/
/*  SET UP VORONOI INPUT FILES                 */
/***********************************************/

/*Generate PBS Script*/
void do_stuff::pbs_full_cats(std::string name, int i){
  std::string test_file;
  std::ofstream res(name.c_str());
  res<<"#!/bin/tcsh\n";
  res<<"#PBS -lwalltime=4:00:00\n";
  res<<"#PBS -lnodes=1:ppn=1\n";
  res<<"#echo $PBS_JOBID\n";
  res<<"cd "+work_dir+"\n";
  res<<"mkdir "+temp_dir+"\n";
  std::ostringstream pos_col_stream,annz_col_stream,lephare_col_stream;
  res<<"if (! -e "+temp_dir+cf_file_name[i]+") cp "+cf_file_path+cf_file_name[i]+" "+temp_dir+"\n";
  if(pos_cols.size()==4){
    for(int k=0;k<4;k++){
      pos_col_stream<<"$"<<pos_cols[k];
      if(k!=3) pos_col_stream<<",";
    }
    res<<"awk 'NR>1 {print "+pos_col_stream.str()+"}' "+temp_dir+cf_file_name[i]+" > "
      +temp_dir+cf_file_name[i]+".incz.cat\n";
    res<<"cp "+temp_dir+cf_file_name[i]+".incz.cat "+cat_dir+"PIPELINE_RUN_"+run+"/NEW_CATALOGUES/\n";
  }
  else{
    for(int k=0;k<3;k++){
      pos_col_stream<<"$"<<pos_cols[k];
      if(k!=2) pos_col_stream<<",";
    }
  }
  if(!strcmp(op_make_annz_test.c_str(),"yes")){
    if(!strcmp(op_spec_z.c_str(),"yes")) annz_col_stream<<"$"<<n_cols+2;
    else annz_col_stream<<"$"<<n_cols+1;
    res<<"if (! -e "+temp_dir+cf_file_name[i]+".annz) cp "+cat_dir+"PIPELINE_RUN_"+run+"/PHOTOZ/"
      +cf_file_name[i]+".annz "+temp_dir+"\n";
    res<<"paste "+temp_dir+cf_file_name[i]+" "+temp_dir+cf_file_name[i]+".annz | awk 'NR>1 {print "
      +pos_col_stream.str()+","+annz_col_stream.str()+"}' > "+temp_dir+cf_file_name[i]+".annz.cat\n";
    res<<"cp "+temp_dir+cf_file_name[i]+".annz.cat "+cat_dir+"PIPELINE_RUN_"+run+"/NEW_CATALOGUES/\n";
  }
  if(!strcmp(op_make_lephare.c_str(),"yes")){
    lephare_col_stream<<"$"<<n_cols+2;
    res<<"if (! -e "+temp_dir+cf_file_name[i]+".lephare) cp "+cat_dir+"PIPELINE_RUN_"+run+"/PHOTOZ/"
      +cf_file_name[i]+".lephare "+temp_dir+"\n";
    res<<"awk 'NR>60' "+temp_dir+cf_file_name[i]+".lephare > "+temp_dir+cf_file_name[i]+".temp\n";
    res<<"paste "+temp_dir+cf_file_name[i]+" "+temp_dir+cf_file_name[i]+".temp | awk 'NR>1 {print "+
      pos_col_stream.str()+","+lephare_col_stream.str()+"}' > "+temp_dir+cf_file_name[i]+".lephare.cat\n";
    res<<"cp "+temp_dir+cf_file_name[i]+".lephare.cat "+cat_dir+"PIPELINE_RUN_"+run+"/NEW_CATALOGUES/\n";
    res<<"rm "+temp_dir+cf_file_name[i]+".temp\n";
  }
  if(!strcmp(rm_temp.c_str(),"yes"))
    res<<"rm "+temp_dir+cf_file_name[i]+" "+temp_dir+cf_file_name[i]+".incz.cat "+temp_dir+cf_file_name[i]+".annz "
      +temp_dir+cf_file_name[i]+".annz.cat "+temp_dir+cf_file_name[i]+".lephare "+temp_dir+cf_file_name[i]+".lephare.cat\n";
  res.close();
}

/*Submit Jobs*/
void do_stuff::make_full_cats(int i,int j){
  char buff[512];
  std::ostringstream idpart,id_stream;
  std::string name1,result,command;
  name1=cf_file_path+cf_file_name[i];
  std::ifstream read_file1(name1.c_str());
  if(read_file1.fail()){
    std::ofstream outlogfile;
    outlogfile.open(log_file_name.c_str(),std::ios_base::app);
    outlogfile<<"- ERROR! [Failed to open file: "<<name1<<"]\n"<<std::flush;
    outlogfile<<"- PIPELINE ABORTED...\n"<<std::flush;
    std::cout<<"- ERROR! [Failed to open file: "<<name1<<"]\n"<<std::flush;
    std::cout<<"- PIPELINE ABORTED...\n"<<std::flush;
    outlogfile.close();
    exit(-1);
  }
  read_file1.close();
  name1=cf_file_path+cf_file_name[i];
  count_cols(name1);
  name1=cat_dir+"PIPELINE_RUN_"+run+"/SCRIPTS/make_full_cats."+cf_file_name[i]+".sh";
  pbs_full_cats(name1,i);
  idpart<<jobid_block[j];
  if(jobid_block[j]==-1) command="qsub "+name1;
  else command="qsub -W depend=afterok:"+idpart.str()+" "+name1;
  FILE *in=popen(command.c_str(),"r");
  fgets(buff, sizeof(buff), in);
  result=buff;
  for(int k=0;k<result.size();k++){
    if(isdigit(result[k]))
      id_stream<<result[k];
  }
  jobid_block[j]=atoi(id_stream.str().c_str());
  pclose(in);
  std::ofstream outlogfile;
  outlogfile.open(log_file_name.c_str(),std::ios_base::app);
  outlogfile<<">> JOB["<<jobid_block[j]<<"]: "<<command<<"\n"<<std::flush;
  outlogfile.close();
  command.clear();
  result.clear();
  idpart.str("");
  id_stream.str("");
}

/***********************************************/
/*  RUN XYZ                                    */
/***********************************************/

/*Generate PBS Script*/
void do_stuff::pbs_xyz(std::string name, int i){
  std::string test_file;
  std::ofstream res(name.c_str());
  res<<"#!/bin/tcsh\n";
  res<<"#PBS -lwalltime=4:00:00\n";
  res<<"#PBS -lnodes=1:ppn=1\n";
  res<<"#echo $PBS_JOBID\n";
  res<<"cd "+work_dir+"\n";
  res<<"mkdir "+temp_dir+"\n";
  if(!strcmp(op_make_full_cats.c_str(),"no"))
    res<<"cp "+cf_file_path+cf_file_name[i]+" "+temp_dir+cf_file_name[i]+"."+op_use_z+".cat\n";
  else
    res<<"if (! -e "+temp_dir+cf_file_name[i]+"."+op_use_z+".cat) cp "+cat_dir+"PIPELINE_RUN_"
      +run+"/NEW_CATALOGUES/"+cf_file_name[i]+"."+op_use_z+".cat "+temp_dir+"\n";
  res<<pipe_path+"run_radec_xyz -i "+temp_dir+cf_file_name[i]+"."+op_use_z+".cat -o "+temp_dir+cf_file_name[i]+".xyz\n";
  res<<"cp "+temp_dir+cf_file_name[i]+".xyz "+cat_dir+"PIPELINE_RUN_"+run+"/VORONOI/\n";
  res.close();
}

/*Submit Jobs*/
void do_stuff::run_xyz(int i,int j){
  char buff[512];
  std::ostringstream idpart,id_stream;
  std::string name,result,command;
  name=cat_dir+"PIPELINE_RUN_"+run+"/SCRIPTS/run_xyz."+cf_file_name[i]+".sh";
  pbs_xyz(name,i);
  idpart<<jobid_block[j];
  if(jobid_block[j]==-1) command="qsub "+name;
  else command="qsub -W depend=afterok:"+idpart.str()+" "+name;
  FILE *in=popen(command.c_str(),"r");
  fgets(buff, sizeof(buff), in);
  result=buff;
  for(int k=0;k<result.size();k++){
    if(isdigit(result[k]))
      id_stream<<result[k];
  }
  jobid_block[j]=atoi(id_stream.str().c_str());
  pclose(in);
  std::ofstream outlogfile;
  outlogfile.open(log_file_name.c_str(),std::ios_base::app);
  outlogfile<<">> JOB["<<jobid_block[j]<<"]: "<<command<<"\n"<<std::flush;
  outlogfile.close();
  command.clear();
  result.clear();
  idpart.str("");
  id_stream.str("");
}

/***********************************************/
/*  FIND XYZ BOX BOUNDARIES                    */
/***********************************************/

/*Generate PBS Script*/
void do_stuff::pbs_find_boundaries(std::string name){
  std::string test_file;
  std::ofstream res(name.c_str());
  res<<"#!/bin/tcsh\n";
  res<<"#PBS -lwalltime=4:00:00\n";
  res<<"#PBS -lnodes=1:ppn=1\n";
  res<<"#echo $PBS_JOBID\n";
  res<<"cd "+work_dir+"\n";
  res<<"mkdir "+temp_dir+"\n";
  for(int o=0;o<cf_file_name.size();o++){
    res<<"if (! -e "+temp_dir+cf_file_name[o]+".xyz) cp "+cat_dir+"PIPELINE_RUN_"+run+"/VORONOI/"
      +cf_file_name[o]+".xyz "+temp_dir+"\n";
  }
  res<<pipe_path+"find_boundaries -i";
  for(int o=0;o<cf_file_name.size();o++){
    res<<" "+temp_dir+cf_file_name[o]+".xyz";
  }
  res<<" -o "+temp_dir+"boundary_file\n";
  res<<"cp "+temp_dir+"boundary_file "+cat_dir+"PIPELINE_RUN_"+run+"/VORONOI/\n";
  res.close();
}

/*Submit Jobs*/
void do_stuff::run_find_boundaries(){
  char buff[512];
  std::ostringstream idpart,id_stream;
  std::string name,result,command;
  name=cat_dir+"PIPELINE_RUN_"+run+"/SCRIPTS/run_find_boundaries.sh";
  pbs_find_boundaries(name);  
  for(int i=0;i<block_size;i++)
    if(jobid_block[i]>-1) idpart<<":"<<jobid_block[i];
  if(jobid_block==neg_block) command="qsub "+name;
  else command="qsub -W depend=afterok"+idpart.str()+" "+name;
  FILE *in=popen(command.c_str(),"r");
  fgets(buff, sizeof(buff), in);
  result=buff;
  for(int k=0;k<result.size();k++){
    if(isdigit(result[k]))
      id_stream<<result[k];
  }
  jobid_block[0]=atoi(id_stream.str().c_str());
  pclose(in);
  std::ofstream outlogfile;
  outlogfile.open(log_file_name.c_str(),std::ios_base::app);
  outlogfile<<">> JOB["<<jobid_block[0]<<"]: "<<command<<"\n"<<std::flush;
  outlogfile.close();
  command.clear();
  result.clear();
  idpart.str("");
  id_stream.str("");
}

/***********************************************/
/*  RUN MAIN SPLIT                             */
/***********************************************/

/*Generate PBS Script*/
void do_stuff::pbs_main_split(std::string name){
  std::string test_file;
  std::ofstream res(name.c_str());
  res<<"#!/bin/tcsh\n";
  res<<"#PBS -lwalltime=4:00:00\n";
  res<<"#PBS -lnodes=1:ppn=1\n";
  res<<"#echo $PBS_JOBID\n";
  res<<"cd "+work_dir+"\n";
  res<<"mkdir "+temp_dir+"\n";
  for(int o=0;o<cf_file_name.size();o++){
    res<<"if (! -e "+temp_dir+cf_file_name[o]+".xyz) cp "+cat_dir+"PIPELINE_RUN_"+run+"/VORONOI/"
      +cf_file_name[o]+".xyz "+temp_dir+"\n";
  }
  res<<"if (! -e "+temp_dir+"boundary_file) cp "+cat_dir+"PIPELINE_RUN_"+run+"/VORONOI/boundary_file "
    +temp_dir+"\n";
  res<<pipe_path+"main_split -i";
  for(int o=0;o<cf_file_name.size();o++){
    res<<" "+temp_dir+cf_file_name[o]+".xyz";
  }
  res<<" -b "+temp_dir+"boundary_file -o "+temp_dir+"boxes_file\n";
  res<<"cp "+temp_dir+"boxes_file "+cat_dir+"PIPELINE_RUN_"+run+"/VORONOI/\n";
  res.close();
}

/*Submit Jobs*/
void do_stuff::run_main_split(){
  char buff[512];
  std::ostringstream idpart,id_stream;
  std::string name,result,command;
  name=cat_dir+"PIPELINE_RUN_"+run+"/SCRIPTS/run_main_split.sh";
  pbs_main_split(name);
  idpart<<jobid_block[0];
  if(jobid_block[0]==-1) command="qsub "+name;
  else command="qsub -W depend=afterok:"+idpart.str()+" "+name;
  FILE *in=popen(command.c_str(),"r");
  fgets(buff, sizeof(buff), in);
  result=buff;
  for(int k=0;k<result.size();k++){
    if(isdigit(result[k]))
      id_stream<<result[k];
  }
  for(int i=0;i<block_size;i++)
    jobid_block[i]=atoi(id_stream.str().c_str());
  pclose(in);
  /*Find N_boxes*/
  int count_sum=0;
  for(int o=0;o<cf_file_name.size();o++){
    std::string name=cf_file_path+cf_file_name[o];
    std::ifstream inFile(name.c_str());
    count_sum+=(int)std::count(std::istreambuf_iterator<char>(inFile),std::istreambuf_iterator<char>(),'\n');  
    inFile.close();
  }
  n_boxes=(int)pow(2,((int)ceil(log((double)count_sum/(double)n_max_split)/log(2.0))));    
  if(n_boxes<=0){
    std::ofstream outlogfile;
    outlogfile.open(log_file_name.c_str(),std::ios_base::app);
    outlogfile<<"ERROR!: MAIN SPLIT ABORT NUMBER OF BOXES <= 0.\n"<<std::flush;
    outlogfile.close();
    exit(-1);
  }
  std::ofstream outlogfile;
  outlogfile.open(log_file_name.c_str(),std::ios_base::app);
  outlogfile<<">> JOB["<<jobid_block[0]<<"]: "<<command<<"\n"<<std::flush;
  outlogfile<<"- NUMBER OF BOXES: "<<n_boxes<<"\n"<<std::flush;
  outlogfile.close();
  command.clear();
  result.clear();
  idpart.str("");
  id_stream.str("");
}

/***********************************************/
/*  RUN SPLIT                                  */
/***********************************************/

/*Generate PBS Script*/
void do_stuff::pbs_split(std::string name,int i,int o){
  std::string test_file;
  std::ofstream res(name.c_str());
  res<<"#!/bin/tcsh\n";
  res<<"#PBS -lwalltime=4:00:00\n";
  res<<"#PBS -lnodes=1:ppn=1\n";
  res<<"#echo $PBS_JOBID\n";
  res<<"cd "+work_dir+"\n";
  res<<"mkdir "+temp_dir+"\n";
  res<<"if (! -e "+temp_dir+cf_file_name[i]+".xyz) cp "+cat_dir+"PIPELINE_RUN_"+run+"/VORONOI/"
    +cf_file_name[i]+".xyz "+temp_dir+"\n";
  res<<"if (! -e "+temp_dir+"boxes_file) cp "+cat_dir+"PIPELINE_RUN_"+run+"/VORONOI/boxes_file "
    +temp_dir+"\n";
  std::ostringstream stream;
  stream<<o;
  res<<pipe_path+"run_split -i "+temp_dir+cf_file_name[i]+".xyz -b "+temp_dir+"boxes_file -n "+stream.str()
    +" -o "+temp_dir+"xyz."+stream.str()+".split";
  res<<"\n";
  res<<"cp "+temp_dir+"xyz."+stream.str()+".split "+cat_dir+"PIPELINE_RUN_"+run+"/VORONOI/SPLITS/\n";
  res.close();
}

/*Submit Jobs*/
void do_stuff::run_split(int i,int j,int o){
  char buff[512];
  std::ostringstream split_stream,idpart,id_stream;
  std::string name,result,command;
  split_stream<<o;
  name=cat_dir+"PIPELINE_RUN_"+run+"/SCRIPTS/SPLITS/run_split."+cf_file_name[i]+"."+split_stream.str()+".sh";
  split_stream.str("");
  pbs_split(name,i,o);
  idpart<<jobid_block[j];
  if(jobid_block[j]==-1) command="qsub "+name;
  else command="qsub -W depend=afterok:"+idpart.str()+" "+name;
  FILE *in=popen(command.c_str(),"r");
  fgets(buff, sizeof(buff), in);
  result=buff;
  for(int k=0;k<result.size();k++){
    if(isdigit(result[k]))
      id_stream<<result[k];
  }
  jobid_block[j]=atoi(id_stream.str().c_str());
  pclose(in);
  std::ofstream outlogfile;
  outlogfile.open(log_file_name.c_str(),std::ios_base::app);
  outlogfile<<">> JOB["<<jobid_block[j]<<"]: "<<command<<"\n"<<std::flush;
  outlogfile.close();
  command.clear();
  result.clear();
  idpart.str("");
  id_stream.str("");
}

/***********************************************/
/*  RUN VORONOI                                */
/***********************************************/

/*Generate PBS Script*/
void do_stuff::pbs_voronoi(std::string name, int i){
  std::string test_file;
  std::ofstream res(name.c_str());
  res<<"#!/bin/tcsh\n";
  res<<"#PBS -lwalltime=4:00:00\n";
  res<<"#PBS -lnodes=1:ppn=1\n";
  res<<"#echo $PBS_JOBID\n";
  res<<"cd "+work_dir+"\n";
  res<<"mkdir "+temp_dir+"\n";
  std::ostringstream stream;
  stream<<i;
  res<<"if (! -e "+temp_dir+"xyz."+stream.str()+".split) cp "+cat_dir+"PIPELINE_RUN_"+run
    +"/VORONOI/SPLITS/xyz."+stream.str()+".split "+temp_dir+"\n";
  res<<"if (! -e "+temp_dir+"boundary_file) cp "+cat_dir+"PIPELINE_RUN_"+run+"/VORONOI/boundary_file "
    +temp_dir+"\n";
  res<<pipe_path+"run_voronoi -i "+temp_dir+"xyz."+stream.str()+".split -b "+temp_dir
    +"boundary_file -o "+temp_dir+"xyz."+stream.str()+".voro\n";
  res<<"cp "+temp_dir+"xyz."+stream.str()+".voro "+cat_dir+"PIPELINE_RUN_"+run+"/VORONOI/SPLITS/\n";
  res.close();
}

void do_stuff::temp_jobs(){
  std::vector<int>::iterator it;
  tempjobs=jobid_block;
  it=std::unique(tempjobs.begin(),tempjobs.end());
  tempjobs.resize(std::distance(tempjobs.begin(),it));
}

/*Submit Jobs*/
void do_stuff::run_voronoi(int i,int j){
  char buff[512];
  std::ostringstream split_stream,idpart,id_stream;
  std::string name,result,command;
  split_stream<<i;
  name=cat_dir+"PIPELINE_RUN_"+run+"/SCRIPTS/SPLITS/run_voronoi."+split_stream.str()+".sh";
  pbs_voronoi(name,i);
  if(i==0){
    temp_jobs();
    for(int o=0;o<tempjobs.size();o++)
      if(tempjobs[o]>0) idpart<<":"<<tempjobs[o];
  }
  else{
    if(j==0) j=1;
    if(j==1) j=0;
    idpart<<":"<<jobid_block[j];
  }
  if(jobid_block==neg_block) command="qsub "+name;
  else command="qsub -W depend=afterok"+idpart.str()+" "+name;
  FILE *in=popen(command.c_str(),"r");
  fgets(buff, sizeof(buff), in);
  result=buff;
  for(int k=0;k<result.size();k++){
    if(isdigit(result[k]))
      id_stream<<result[k];
  }
  jobid_block[j]=atoi(id_stream.str().c_str());
  pclose(in);
  std::ofstream outlogfile;
  outlogfile.open(log_file_name.c_str(),std::ios_base::app);
  outlogfile<<">> JOB["<<jobid_block[j]<<"]: "<<command<<"\n"<<std::flush;
  outlogfile.close();
  command.clear();
  result.clear();
  idpart.str("");
  id_stream.str("");
}

/***********************************************/
/*  RUN CLEAN UP                               */
/***********************************************/

/*Generate PBS Script*/
void do_stuff::pbs_cleanup(std::string name,int i,int k){
  std::string test_file;
  std::ofstream res(name.c_str());
  res<<"#!/bin/tcsh\n";
  res<<"#PBS -lwalltime=4:00:00\n";
  res<<"#PBS -lnodes=1:ppn=1\n";
  res<<"#echo $PBS_JOBID\n";
  res<<"cd "+work_dir+"\n";
  res<<"mkdir "+temp_dir+"\n";
  std::ostringstream stream1,stream2;
  stream1<<i;
  stream2<<k;
  res<<"if (! -e "+temp_dir+"xyz."+stream1.str()+".voro) cp "+cat_dir+"PIPELINE_RUN_"+run
    +"/VORONOI/SPLITS/xyz."+stream1.str()+".voro "+temp_dir+"\n";
  res<<"if (! -e "+temp_dir+"xyz."+stream2.str()+".voro) cp "+cat_dir+"PIPELINE_RUN_"+run
    +"/VORONOI/SPLITS/xyz."+stream2.str()+".voro "+temp_dir+"\n";
  res<<pipe_path+"run_cleanup -i "+temp_dir+"xyz."+stream1.str()+".voro "+temp_dir+"xyz."+stream2.str()+".voro\n";
  res<<"cp "+temp_dir+"xyz."+stream1.str()+".voro "+temp_dir+"xyz."+stream2.str()+".voro "+cat_dir+"PIPELINE_RUN_"+run
    +"/VORONOI/SPLITS/\n";
  res.close();
}

/*Submit Jobs*/
void do_stuff::run_cleanup(int i,int j,int k){
  char buff[512];
  std::ostringstream split_stream1,split_stream2,idpart,id_stream;
  std::string name,result,command;
  split_stream1<<i;
  split_stream2<<k;
  name=cat_dir+"PIPELINE_RUN_"+run+"/SCRIPTS/SPLITS/run_cleanup."+split_stream1.str()+"_"+split_stream2.str()+".sh";
  pbs_cleanup(name,i,k);
  if(i==0 && j==0 && k==1) temp_jobs();      
  if(j==(k-1)){
    for(int o=0;o<tempjobs.size();o++)
      if(tempjobs[o]>0) idpart<<":"<<tempjobs[o];
  }
  else
    idpart<<":"<<jobid_block[j];
  if(jobid_block==neg_block) command="qsub "+name;
  else command="qsub -W depend=afterok"+idpart.str()+" "+name;
  FILE *in=popen(command.c_str(),"r");
  fgets(buff, sizeof(buff), in);
  result=buff;
  for(int kk=0;kk<result.size();kk++){
    if(isdigit(result[kk]))
      id_stream<<result[kk];
  }
  jobid_block[j]=atoi(id_stream.str().c_str());
  pclose(in);
  std::ofstream outlogfile;
  outlogfile.open(log_file_name.c_str(),std::ios_base::app);
  outlogfile<<">> JOB["<<jobid_block[j]<<"]: "<<command<<"\n"<<std::flush;
  outlogfile.close();
  command.clear();
  result.clear();
  idpart.str("");
  id_stream.str("");
}

/***********************************************/
/*  RUN PIPELINKING LINKING                    */
/***********************************************/

/*Generate PBS Script*/
void do_stuff::pbs_pipeline_linking(std::string name, int i){
  std::string test_file;
  std::ofstream res(name.c_str());
  res<<"#!/bin/tcsh\n";
  res<<"#PBS -lwalltime=4:00:00\n";
  res<<"#PBS -lnodes=1:ppn=1\n";
  res<<"#echo $PBS_JOBID\n";
  res<<"cd "+work_dir+"\n";
  res<<"mkdir "+temp_dir+"\n";
  std::ostringstream stream;
  stream<<i;
  res<<"if (! -e "+temp_dir+"xyz."+stream.str()+".voro) cp "+cat_dir+"PIPELINE_RUN_"+run
    +"/VORONOI/SPLITS/xyz."+stream.str()+".voro "+temp_dir+"\n";
  res<<pipe_path+"run_pipeline_linking -i "+temp_dir+"xyz."+stream.str()+".voro -f "+cf_file_list+" -pre "
    +temp_dir+" -post ."+op_use_z+".cat -o "+temp_dir+"xyz."+stream.str()+".link\n";
  res<<"cp "+temp_dir+"xyz."+stream.str()+".link "+cat_dir+"PIPELINE_RUN_"+run+"/VORONOI/SPLITS/\n";
  res.close();
}

/*Submit Jobs*/
void do_stuff::run_pipeline_linking(int i,int j){
  char buff[512];
  std::ostringstream split_stream,idpart,id_stream;
  std::string name,result,command;
  split_stream<<i;
  name=cat_dir+"PIPELINE_RUN_"+run+"/SCRIPTS/SPLITS/run_pipeline_linking."+split_stream.str()+".sh";
  pbs_pipeline_linking(name,i);
  if(j==0) temp_jobs();
  if(i==j){
    for(int o=0;o<tempjobs.size();o++)
      if(tempjobs[o]>0) idpart<<":"<<tempjobs[o];
  }
  else
    idpart<<":"<<jobid_block[j];
  if(jobid_block==neg_block) command="qsub "+name;
  else command="qsub -W depend=afterok"+idpart.str()+" "+name;
  FILE *in=popen(command.c_str(),"r");
  fgets(buff, sizeof(buff), in);
  result=buff;
  for(int k=0;k<result.size();k++){
    if(isdigit(result[k]))
      id_stream<<result[k];
  }
  jobid_block[j]=atoi(id_stream.str().c_str());
  pclose(in);
  std::ofstream outlogfile;
  outlogfile.open(log_file_name.c_str(),std::ios_base::app);
  outlogfile<<">> JOB["<<jobid_block[j]<<"]: "<<command<<"\n"<<std::flush;
  outlogfile.close();
  command.clear();
  result.clear();
  idpart.str("");
  id_stream.str("");
}

/***********************************************/
/*  RUN DENSITY CUTS                           */
/***********************************************/

/*Generate PBS Script*/
void do_stuff::pbs_density_cuts(std::string name, int i){
  std::string test_file;
  std::ofstream res(name.c_str());
  res<<"#!/bin/tcsh\n";
  res<<"#PBS -lwalltime=4:00:00\n";
  res<<"#PBS -lnodes=1:ppn=1\n";
  res<<"#echo $PBS_JOBID\n";
  res<<"cd "+work_dir+"\n";
  res<<"mkdir "+temp_dir+"\n";
  std::ostringstream stream;
  stream<<i;
  res<<"if (! -e "+temp_dir+"xyz."+stream.str()+".link) cp "+cat_dir+"PIPELINE_RUN_"+run
    +"/VORONOI/SPLITS/xyz."+stream.str()+".link "+temp_dir+"\n";
  res<<pipe_path+"density_cuts -i "+temp_dir+"xyz."+stream.str()+".link -c "+density_cut_value+" -o "
    +temp_dir+"xyz."+stream.str()+".cut\n";
  res<<"cp "+temp_dir+"xyz."+stream.str()+".cut "+cat_dir+"PIPELINE_RUN_"+run+"/VORONOI/SPLITS/\n";
  res.close();
}

/*Submit Jobs*/
void do_stuff::run_density_cuts(int i,int j){
  char buff[512];
  std::ostringstream split_stream,idpart,id_stream;
  std::string name,result,command;
  split_stream<<i;
  name=cat_dir+"PIPELINE_RUN_"+run+"/SCRIPTS/SPLITS/run_density_cuts."+split_stream.str()+".sh";
  pbs_density_cuts(name,i);  
  idpart<<jobid_block[j];
  if(jobid_block[j]==-1) command="qsub "+name;
  else command="qsub -W depend=afterok:"+idpart.str()+" "+name;
  FILE *in=popen(command.c_str(),"r");
  fgets(buff, sizeof(buff), in);
  result=buff;
  for(int k=0;k<result.size();k++){
    if(isdigit(result[k]))
      id_stream<<result[k];
  }
  //jobid_block[j]=atoi(id_stream.str().c_str());
  for(int i=0;i<block_size;i++)
    jobid_block[i]=atoi(id_stream.str().c_str());
  pclose(in);
  std::ofstream outlogfile;
  outlogfile.open(log_file_name.c_str(),std::ios_base::app);
  outlogfile<<">> JOB["<<jobid_block[j]<<"]: "<<command<<"\n"<<std::flush;
  outlogfile.close();
  command.clear();
  result.clear();
  idpart.str("");
  id_stream.str("");
}

/***********************************************/
/*  RUN SPLIT RADEC                            */
/***********************************************/

/*Generate PBS Script*/
void do_stuff::pbs_split_radec(std::string name,int i,int o){
  std::string test_file;
  std::ofstream res(name.c_str());
  res<<"#!/bin/tcsh\n";
  res<<"#PBS -lwalltime=4:00:00\n";
  res<<"#PBS -lnodes=1:ppn=1\n";
  res<<"#echo $PBS_JOBID\n";
  res<<"cd "+work_dir+"\n";
  res<<"mkdir "+temp_dir+"\n";
  res<<"source /home/cosmos/sourceBinPATHS\n";
  res<<"source /home/cosmos/sourceLibPATHS\n";
  res<<"source /home/cosmos/sourceIntelCompilers\n";
  std::ostringstream stream1,stream2,pix,perc;
  stream1<<i;
  stream2<<o;
  pix<<num_pixels;
  perc<<inter_perc;
  res<<"cp "+cat_dir+"PIPELINE_RUN_"+run+"/VORONOI/SPLITS/xyz."+stream1.str()+".cut "+temp_dir+"\n";
  res<<pipe_path+"run_splitradec2 -i "+temp_dir+"xyz."+stream1.str()+".cut -o "+temp_dir+"radec."+stream2.str()+".split -p "
    +pix.str()+" -inters "+perc.str()+" -n "+stream2.str()+"\n";
  res<<"cp "+temp_dir+"radec."+stream2.str()+".split "+cat_dir+"PIPELINE_RUN_"+run+"/VORONOI/SPLITS_RADEC/\n";
  res.close();
}

///  res<<"find "+temp_dir+" -size 0 -delete\n";

/*Submit Jobs*/
void do_stuff::run_split_radec(int i,int j,int o){
  char buff[512];
  std::ostringstream split_stream1,split_stream2,idpart,id_stream;
  std::string name,result,command;
  split_stream1<<i;
  split_stream2<<o;
  name=cat_dir+"PIPELINE_RUN_"+run+"/SCRIPTS/SPLITS_RADEC/run_split_radec."+split_stream1.str()+"."+split_stream2.str()+".sh";
  split_stream1.str("");
  split_stream2.str("");
  pbs_split_radec(name,i,o);
  idpart<<jobid_block[j];
  if(jobid_block[j]==-1) command="qsub "+name;
  else command="qsub -W depend=afterok:"+idpart.str()+" "+name;
  FILE *in=popen(command.c_str(),"r");
  fgets(buff, sizeof(buff), in);
  result=buff;
  for(int k=0;k<result.size();k++){
    if(isdigit(result[k]))
      id_stream<<result[k];
  }
  jobid_block[j]=atoi(id_stream.str().c_str());
  pclose(in);
  std::ofstream outlogfile;
  outlogfile.open(log_file_name.c_str(),std::ios_base::app);
  outlogfile<<">> JOB["<<jobid_block[j]<<"]: "<<command<<"\n"<<std::flush;
  outlogfile.close();
  command.clear();
  result.clear();
  idpart.str("");
  id_stream.str("");
}

/***********************************************/

int factorial(int num);

/***MAIN***/

int main(int argc, char *argv[]){

  /***********************************************/

  if(argc!=2){
    std::cout<<"Please provide parameters file (e.g. params.ini) path as argument 1.\n"<<std::flush;
    exit(-1);
  }
  class do_stuff do_it;
  int start,end,job_pos,counter;

  /***********************************************/
  /*  READ PARAMETER FILE                        */
  /***********************************************/

  do_it.read_params(argv[1]);
  
  /***********************************************/
  /*  SET UP RUN DIRECTORIES AND PIPELINE LOG    */
  /***********************************************/

  if(!strcmp(do_it.op_dir_set_up.c_str(),"yes")){
    do_it.set_up(do_it.run_number);
    do_it.dir_set_up();
    do_it.set_up_log();
  }

  /***********************************************/
  /*  READ LIST OF CATALOGUE FILES               */
  /***********************************************/

  do_it.read_file_list();

  /***********************************************/
  /*  SET UP ANNZ AND LEPHARE INPUT FILES        */
  /***********************************************/

  if(!strcmp(do_it.op_make_initial_cats.c_str(),"yes")){
    counter=0;
    for(;;){ /*Loop over original files*/
      start=counter;
      end=counter+do_it.block_size;
      int j=0;
      for(int i=start;i<end;i++){
	do_it.make_initial_cats(i,j);
	counter++; j++;
	if(counter>=do_it.cf_file_name.size()) break;
      }
      if(counter>=do_it.cf_file_name.size()) break;
    }
  }

  /***********************************************/
  /*  TRAIN ANNZ NETWORK                         */
  /***********************************************/
  
  if(!strcmp(do_it.op_make_annz_train.c_str(),"yes")){
    counter=0;
    for(;;){ /*Loop over number of networks*/
      start=counter;
      end=counter+do_it.block_size;
      int j=0;
      for(int i=start;i<end;i++){
	do_it.make_annz_train(i,j);
	counter++; j++;
	if(counter>=do_it.num_nets) break;
      }
      if(counter>=do_it.num_nets) break;
    }
  }

  /***********************************************/
  /*  RUN ANNZ                                   */
  /***********************************************/
  
  if(!strcmp(do_it.op_make_annz_test.c_str(),"yes")){
    counter=0;
    for(;;){ /*Loop over original files*/
      start=counter;
      end=counter+do_it.block_size;
      int j=0;
      for(int i=start;i<end;i++){
	do_it.make_annz_test(i,j);
	counter++; j++;
	if(counter>=do_it.cf_file_name.size()) break;
      }
      if(counter>=do_it.cf_file_name.size()) break;
    }
  }
  
  /***********************************************/
  /*  RUN LEPHARE                                */
  /***********************************************/

  if(!strcmp(do_it.op_make_lephare.c_str(),"yes")){
    counter=0;
    for(;;){ /*Loop over original files*/
      start=counter;
      end=counter+do_it.block_size;
      int j=0;
      for(int i=start;i<end;i++){
	do_it.make_lephare(i,j);
	counter++; j++;
	if(counter>=do_it.cf_file_name.size()) break;
      }
      if(counter>=do_it.cf_file_name.size()) break;
    }
  }

  /***********************************************/
  /*  SET UP VORONOI INPUT FILES                 */
  /***********************************************/

  if(!strcmp(do_it.op_make_full_cats.c_str(),"yes")){
    counter=0;
    for(;;){ /*Loop over original files*/
      start=counter;
      end=counter+do_it.block_size;
      int j=0;
      for(int i=start;i<end;i++){
	do_it.make_full_cats(i,j);
	counter++; j++;
	if(counter>=do_it.cf_file_name.size()) break;
      }
      if(counter>=do_it.cf_file_name.size()) break;
    }
  }

  /***********************************************/
  /*  RUN XYZ                                    */
  /***********************************************/
  
  if(!strcmp(do_it.op_run_xyz.c_str(),"yes")){
    counter=0;
    for(;;){ /*Loop over original files*/
      start=counter;
      end=counter+do_it.block_size;
      int j=0;
      for(int i=start;i<end;i++){
	do_it.run_xyz(i,j);
	counter++; j++;
	if(counter>=do_it.cf_file_name.size()) break;
      }
      if(counter>=do_it.cf_file_name.size()) break;
    }
  }

  /***********************************************/
  /*  FIND XYZ BOX BOUNDARIES                    */
  /***********************************************/

  if(!strcmp(do_it.op_run_find_boundaries.c_str(),"yes")) do_it.run_find_boundaries();
  
  /***********************************************/
  /*  RUN MAIN SPLIT                             */
  /***********************************************/

  if(!strcmp(do_it.op_run_main_split.c_str(),"yes")) do_it.run_main_split();

  /***********************************************/
  /*  RUN SPLIT                                  */
  /***********************************************/

  if(!strcmp(do_it.op_run_split.c_str(),"yes")){
    counter=0; 
    int i=0,o=0;
    for(;;){ /*Loop over original files and split boxes*/
      start=counter;
      end=counter+do_it.block_size;
      int j=0;
      for(int k=start;k<end;k++){
	do_it.run_split(i,j,o);
	j++; o++;
	if(o==do_it.n_boxes){
	  i++; o=0; counter++;
	}
	if(counter>=do_it.cf_file_name.size()) break;
      }
      if(counter>=do_it.cf_file_name.size()) break;
    }  
  }

  /***********************************************/
  /*  RUN VORONOI                                */
  /***********************************************/

  if(!strcmp(do_it.op_run_voronoi.c_str(),"yes")){
    counter=0;
    for(;;){ /*Loop over split files*/
      start=counter;
      end=counter+do_it.block_size;
      int j=0;
      for(int i=start;i<end;i++){
	do_it.run_voronoi(i,j);
	counter++; j++;
	if(counter>=do_it.n_boxes) break;
      }
      if(counter>=do_it.n_boxes) break;
    }
  }

  /***********************************************/
  /*  RUN CLEAN UP                               */
  /***********************************************/

  if(!strcmp(do_it.op_run_cleaup.c_str(),"yes")){
    int i_count=0; counter=1;
    for(;;){ /*Loop over split files*/
      start=counter;
      end=counter+do_it.block_size;
      int j=0;
      for(int k=start;k<end;k++){
	do_it.run_cleanup(i_count,j,k);
	counter++;
	j++;
	if(i_count>=do_it.n_boxes-1) break;
	if(counter>=do_it.n_boxes){
	  i_count++;
	  counter=i_count+1;
	  break;
	}
      }
      if(i_count>=do_it.n_boxes-1) break;
    }
  }

  /***********************************************/
  /*  RUN PIPELINE LINKING                       */
  /***********************************************/
 
  if(!strcmp(do_it.op_run_pipeline_linking.c_str(),"yes")){
    counter=0;
    for(;;){ /*Loop over split files*/
      start=counter;
      end=counter+do_it.block_size;
      int j=0;
      for(int i=start;i<end;i++){
	do_it.run_pipeline_linking(i,j);
	counter++; j++;
	if(counter>=do_it.n_boxes) break;
      }
      if(counter>=do_it.n_boxes) break;
    }
  }

  /***********************************************/
  /*  RUN DENSITY CUTS                           */
  /***********************************************/
  
  if(!strcmp(do_it.op_density_cuts.c_str(),"yes")){
    counter=0;
    for(;;){ /*Loop over split files*/
      start=counter;
      end=counter+do_it.block_size;
      int j=0;
      for(int i=start;i<end;i++){
	do_it.run_density_cuts(i,j);
	counter++; j++;
	if(counter>=do_it.n_boxes) break;
      }
      if(counter>=do_it.n_boxes) break;
    }
  }

  /***********************************************/
  /*  RUN SPLIT RADEC                            */
  /***********************************************/

  if(!strcmp(do_it.op_run_split_radec.c_str(),"yes")){
    counter=0; 
    int i=0,o=0;
    for(;;){ /*Loop over original files and split boxes*/
      start=counter;
      end=counter+do_it.block_size;
      int j=0;
      for(int k=start;k<end;k++){
	do_it.run_split_radec(i,j,o);
	j++; o++;      
	if(o==do_it.num_pixels){
	  i++; o=0; counter++;
	}
	if(counter>=do_it.n_boxes) break;
      }
      if(counter>=do_it.n_boxes) break;
    }  
  }
  
  /***********************************************/

  return 0;
}

/***END MAIN***/

/***********************************************/
/*  FACTORIAL FUNCTION                         */
/***********************************************/

int factorial(int num){
  int result=1;
  for (int i=1;i<=num;++i) result=result*=i;
  return result;
}

/***********************************************/
