char buff[512];
string script_name, command, result;
std::ostringstream id_stream;
int job_id;                                 

script_name = "job_script.pbs"                   //Name of PBS script to be submitted.

command = "qsub "+script_name;                   //Define job submission command line.
FILE *in = popen(command.c_str(), "r");          //Execute command line.
fgets(buff, sizeof(buff), in);                   //Retrieve job ID from stdout and store in buff character.
result = buff;                                   //Store character contents in result string.
for(int i = 0; i < result.size(); i++)           //Loop through string elements.
  if(isdigit(result[i])) id_stream << result[i]; //If string element is a digit, store in id_stream.
job_id = atoi(id_stream.str().c_str());          //Store job ID as an integer.
