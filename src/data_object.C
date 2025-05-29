//data_object.C
//Original Author: Ezekiel Wertz
//Update Author: Kate Evans
//Companion implementation. Class to represent a data object file for GMn/nTPE analysis. Anything we need to know about data files.

#include "../include/data_object.h"
#include <iostream>
#include "TMath.h"

//private helper function. Requires that private class variables are initialized first
//need to test this function carefully. Should be correct for data file directory structure
TString data_object::makeInputFileName(){
  string input_directory = "/volatile/halla/sbs/sbs-gen/GEN_REPLAYS/pass2/TEST/try8";
  TString inputfile;
  const char *input_directory_char = input_directory.c_str();
  const char *pass_char = pass.Data();
  const char *kin_char = kinematic.Data();
  const char *tar_char = target.Data();
  if(pass == "pass2"){

    if(kinematic == "GEN2"){
      if(target == "He3"){
	inputfile = Form("%s/%s/rootfiles/e1209016_fullreplay_%i_*.root",input_directory_char,kin_char,run);
      } // end if He3
      if(target == "H2"){
	inputfile = Form("%s/%s/%s/rootfiles/e1209016_fullreplay_%i_*.root",input_directory_char,kin_char,tar_char,run);
      } // end if H2
    } // end if GEN2

    if(kinematic == "GEN3"){
      if(target == "He3"){
	inputfile = Form("%s/%s/rootfiles/e1209016_fullreplay_%i_*.root",input_directory_char,kin_char,run);
      } // end if He3
      if(target == "H2"){
	inputfile = Form("%s/%s/%s/rootfiles/e1209016_fullreplay_%i_*.root",input_directory_char,kin_char,tar_char,run);
      } // end if H2
    } // end if GEN3

    if(kinematic == "GEN4a"){
      if(target == "He3"){
	inputfile = Form("%s/%s/rootfiles/e1209016_fullreplay_%i_*.root",input_directory_char,kin_char,run);
      } // end if He3
    } // end if GEN4a

    if(kinematic == "GEN4b"){
      if(target == "He3"){
	inputfile = Form("%s/%s/rootfiles/e1209016_fullreplay_%i_*.root",input_directory_char,kin_char,run);
      } // end if He3
    } // end if GEN4b

  } // end if pass2
  else{
    //make an error. Some how we got not pass 0,1, or 2
    cout << "Error: Pass variable was given that is not pass 0,1, or 2 " << pass << " !" << endl;
  }
return inputfile;
} // end makeInputFileName

data_object::data_object(int runnum, const char *data_file_name, const char *kinematic_file_name, TString Kin, TString SBS_field, TString targ, TString daPass){

ifstream datafile(data_file_name);

if(datafile.fail()){
  cout << "Error:There was a problem with the data file " << data_file_name << ". Figure it out nerd!" << endl;
  return;
}

TString currentLine;
bool gotRun = false;
TString runnum_string = utility::intToTString(runnum);

while(currentLine.ReadLine(datafile)){
  if(currentLine.BeginsWith(runnum_string)){
    TObjArray *tokens = currentLine.Tokensize("  ");
    run = (((TObjString*) (*tokens)[0])->GetString()).Atoi();
    pass = ((TObjString*) (*tokens)[1])->GetString();
    kinematic = ((TObjString*) (*tokens)[2])->GetString();
    target = ((TObjString*) (*tokens)[3])->GetString();
    sbs_field = (((TObjString*) (*tokens)[4])->GetString()).Atoi();
    // define other things from the data map?
  } // end if line begins with runnum
} // end while reading lines

} // end data_object
