//****************************************//
//***            CATFORM               ***//
//****************************************//
//*** Written by Samuel Farrens (2012) ***//
//*** LAST UPDATE: 20-11-2012          ***//
//*** Contact: farrens@ieec.uab.es     ***//
//****************************************//

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
#include<sys/stat.h>
#include <unistd.h>

class do_stuff{
private:
public:
  int block_size;
  std::string user,work_dir,cat_dir,pbs_dir,temp_dir,tab_path,node_list;
  std::vector<int> jobid_block;
  std::vector<std::string> nodes,file_type,og_file_list,og_file_path;
  std::vector<std::vector<std::string> > og_file_name,columns;
  void read_params(char *file_named);
  void read_file_list();
  void split(const std::string &str,std::vector<std::string> &tokens,const std::string &delimiter);
  void whileEOF(std::ifstream &inputfile,std::vector<std::string> &header,std::vector<std::string> &values,
		const std::string &comment_str);
  void set_up();
  void pbs_initial_cats(std::string name, int i);
  void make_initial_cats(int i,int j);
  void pbs_remove_temp_files(std::string name);
  void remove_temp_files();
};

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

void do_stuff::read_params(char *file_named){
  int colcount=-1;
  std::ifstream file(file_named);
  std::vector<std::string> header,lines,values;
  whileEOF(file,header,lines,"#");
  file.close();
  for(int i=0;i<lines.size();++i){
    values.clear();
    split(lines[i],values," "); 
    if(values[0]=="user") user=values[1];
    if(values[0]=="work_dir") work_dir=values[1];
    if(values[0]=="cat_dir") cat_dir=values[1];
    if(values[0]=="pbs_dir") pbs_dir=values[1];
    if(values[0]=="temp_dir") temp_dir=values[1];
    if(values[0]=="tab_path") tab_path=values[1];
    if(values[0]=="block_size") block_size=atoi(values[1].c_str());
    if(values[0]=="file_type") file_type.push_back(values[1]);
    if(values[0]=="node_list") node_list=values[1];
    if(values[0]=="og_file_path") og_file_path.push_back(values[1]);
    if(values[0]=="og_file_list"){
      og_file_list.push_back(values[1]);
      if(values.size()>2){
	columns.push_back(std::vector<std::string>());
	colcount++;
      	for(int j=2;j<values.size();j++){
	  columns[colcount].push_back(values[j]);
	}
      }
    }
  }
  std::ifstream test_file(node_list.c_str());
  if(test_file.fail()){
    std::cout<<"ERROR! "<<node_list<<" NOT FOUND.\n"<<std::flush;
    std::cout<<"ABORT...\n"<<std::flush;
    exit(-1);
  }
}

void do_stuff::read_file_list(){
  std::string s_var;
  for(int i=0;i<og_file_list.size();i++){
    og_file_name.push_back(std::vector<std::string>());
    std::ifstream read_file(og_file_list[i].c_str()); 
    for(;;){
      read_file>>s_var;
      if(read_file.eof()) break;
      og_file_name[i].push_back(s_var);
    }
    read_file.close();
  }
}

void do_stuff::set_up(){
  for(int i=0;i<block_size;i++) jobid_block.push_back(-1);
  if(temp_dir.empty()) temp_dir=user+"_CATFORM_TEMPORARY_FILES/";
  mkdir(cat_dir.c_str(),S_IRWXU|S_IRGRP|S_IXGRP);
  mkdir(pbs_dir.c_str(),S_IRWXU|S_IRGRP|S_IXGRP);
}

void do_stuff::pbs_initial_cats(std::string name, int i){
  std::ofstream res(name.c_str());
  res<<"#!/bin/tcsh\n";
  res<<"#PBS -lwalltime=4:00:00\n";
  res<<"#PBS -lnodes=1:ppn=1\n";
  res<<"#echo $PBS_JOBID\n";
  res<<"cd "+work_dir+"\n";
  res<<"mkdir "+temp_dir+"\n";
  std::ostringstream stream,col_stream,paste_stream;
  stream<<i;
  for(int o=0;o<og_file_list.size();o++){
    res<<"cp "+og_file_path[o]+og_file_name[o][i]+" "+temp_dir+"\n";
    paste_stream<<temp_dir<<og_file_name[o][i]<<".cat ";
    for(int k=0;k<columns[o].size();k++){
      if(!strcmp(file_type[o].c_str(),"fits")) col_stream<<columns[o][k];
      else if(!strcmp(file_type[o].c_str(),"ascii")) col_stream<<"$"<<columns[o][k];
      if(k!=columns[o].size()-1){
	if(!strcmp(file_type[o].c_str(),"fits")) col_stream<<";";
	else if(!strcmp(file_type[o].c_str(),"ascii")) col_stream<<",";
      }
    }
    if(!strcmp(file_type[o].c_str(),"fits"))
      res<<tab_path+"tablist '"+temp_dir+og_file_name[o][i]+"[col "+col_stream.str()+"]' > "+temp_dir+og_file_name[o][i]+".cat\n";
    else if(!strcmp(file_type[o].c_str(),"ascii"))
      res<<"awk '{print "+col_stream.str()+"}' "+temp_dir+og_file_name[o][i]+" > "+temp_dir+og_file_name[o][i]+".cat\n";
    col_stream.str("");
  }
  if(og_file_list.size()>1){
    res<<"paste "+paste_stream.str()+"> "+temp_dir+og_file_name[0][i]+".total.cat\n";
    res<<"cp "+temp_dir+og_file_name[0][i]+".total.cat "+cat_dir+"\n";
  }
  else
    res<<"cp "+temp_dir+og_file_name[0][i]+".cat "+cat_dir+"\n";
  res.close();
}

void do_stuff::make_initial_cats(int i,int j){
  char buff[512];
  std::ostringstream id_stream,idpart;
  std::string name,name2,command,result;
  name=cat_dir+og_file_name[0][i]+".cat";
  name2=cat_dir+og_file_name[0][i]+".total.cat";
  std::ifstream read_file(name.c_str());
  std::ifstream read_file2(name2.c_str());
  if(read_file.fail() && read_file2.fail()){
    name=pbs_dir+og_file_name[0][i]+".sh";
    pbs_initial_cats(name,i);
    idpart<<jobid_block[j];
    if(jobid_block[j]==-1) command="qsub "+name;
    else command="qsub -W depend=afterok:"+idpart.str()+" "+name;
    std::cout<<command<<"\n"<<std::flush;
    FILE *in=popen(command.c_str(),"r");
    fgets(buff, sizeof(buff), in);
    result=buff;
    for(int k=0;k<result.size();k++){
      if(isdigit(result[k]))
	id_stream<<result[k];
    }
    jobid_block[j]=atoi(id_stream.str().c_str());
    pclose(in);
  }
  read_file.close();
  read_file2.close();
}

void do_stuff::pbs_remove_temp_files(std::string name){
  std::ofstream res(name.c_str());
  res<<"#!/bin/tcsh\n";
  res<<"#PBS -lwalltime=4:00:00\n";
  res<<"#PBS -lnodes=1:ppn=1\n";
  res<<"#echo $PBS_JOBID\n";
  res<<"cd "+work_dir+"\n";
  res<<"rm -r "+temp_dir+"\n";
  res.close();
}

void do_stuff::remove_temp_files(){
  std::string s_var,name,command;
  std::ifstream read_node(node_list.c_str());
  for(;;){
    read_node>>s_var;
    if(read_node.eof()) break;
    nodes.push_back(s_var);
  }
  read_node.close();
  std::ostringstream idpart;
  name=pbs_dir+"remove_temp_files.sh";
  for(int i=0;i<block_size;i++)
    if(jobid_block[i]>-1) idpart<<":"<<jobid_block[i];
  if(!idpart.str().empty()){
    pbs_remove_temp_files(name);
    for(int j=0;j<nodes.size();j++){
      command="qsub -l nodes="+nodes[j]+" -W depend=afterok"+idpart.str()+" "+name;
      std::cout<<command<<"\n"<<std::flush;
      FILE *in=popen(command.c_str(),"r");
      pclose(in);
    }
  }
}

int main(int argc, char *argv[]){
  if(argc!=2){
    std::cout<<"Please provide parameters file (e.g. cf_params.ini) path as argument 1.\n"<<std::flush;
    exit(-1);
  }
  class do_stuff do_it;
  do_it.read_params(argv[1]);
  do_it.set_up();
  do_it.read_file_list();  
  int start,end,counter=0;
  for(;;){ /*Loop over original files*/
    start=counter;
    end=counter+do_it.block_size;
    int j=0;
    for(int i=start;i<end;i++){
      do_it.make_initial_cats(i,j);
      counter++; j++;
      if(counter>=do_it.og_file_name[0].size()) break;
    }
    if(counter>=do_it.og_file_name[0].size()) break;
  }
  do_it.remove_temp_files();
  return 0;
}
