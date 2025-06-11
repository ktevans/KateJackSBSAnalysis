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
  inputfile = "";
/*
  const char *input_directory_char = input_directory.c_str();
  const char *pass_char = pass.Data();
  const char *kin_char = kinematic.Data();
  const char *tar_char = target.Data();
  if(pass == "pass2"){

    if(kinematic == "GEN2"){
      if(target == "He3"){
	inputfile = Form("%s/%s/rootfiles/e1209016_fullreplay_%i_*.root",input_directory_char,kin_char,run);
      } // end if He3
      else if(target == "H2"){
	inputfile = Form("%s/%s/%s/rootfiles/e1209016_fullreplay_%i_*.root",input_directory_char,kin_char,tar_char,run);
      } // end if H2
      else {
	cout << "ERROR!! There is no target: " << target << ", dummy!" << endl;
      }
    } // end if GEN2

    else if(kinematic == "GEN3"){
      if(target == "He3"){
	inputfile = Form("%s/%s/rootfiles/e1209016_fullreplay_%i_*.root",input_directory_char,kin_char,run);
      } // end if He3
      else if(target == "H2"){
	inputfile = Form("%s/%s/%s/rootfiles/e1209016_fullreplay_%i_*.root",input_directory_char,kin_char,tar_char,run);
      } // end if H2
      else {
	cout << "ERROR!! There is no target: " << target << ", dummy!" << endl;
      }
    } // end if GEN3

    else if(kinematic == "GEN4a"){
      if(target == "He3"){
	inputfile = Form("%s/%s/rootfiles/e1209016_fullreplay_%i_*.root",input_directory_char,kin_char,run);
      } // end if He3
      else {
	cout << "ERROR!! There is no target: " << target << ", dummy!" << endl;
      }
    } // end if GEN4a

    else if(kinematic == "GEN4b"){
      if(target == "He3"){
	inputfile = Form("%s/%s/rootfiles/e1209016_fullreplay_%i_*.root",input_directory_char,kin_char,run);
      } // end if He3
      else {
	cout << "ERROR!! There is no target: " << target << ", dummy!" << endl;
      }
    } // end if GEN4b

    else {
      cout << "ERROR!! There is no kinematic: " << kinematic << ", dummy!" << endl;
    }

  } // end if pass2
  else{
    //make an error. Some how we got not pass 0,1, or 2
    cout << "Error: Pass variable was given that is not pass 0,1, or 2. It's actually pass " << pass << "! That's not a pass. Stupid." << endl;
  } */
return inputfile;
} // end makeInputFileName

data_object::data_object(int runnum, const char *data_file_name, const char *kinematic_file_name, TString Kin, TString SBS_field, TString targ, TString daPass){
/*
ifstream datafile(data_file_name);

if(datafile.fail()){
  cout << "Error:There was a problem with the data file " << data_file_name << ". Figure it out, nerd!" << endl;
  return;
}

TString currentLine;
bool gotRun = false;
TString runnum_string = utility::intToTString(runnum);

while(currentLine.ReadLine(datafile)){
  if(currentLine.BeginsWith(runnum_string)){
    TObjArray *tokens = currentLine.Tokenize(",");
    run = (((TObjString*) (*tokens)[0])->GetString()).Atoi();
    pass = ((TObjString*) (*tokens)[1])->GetString();
    kinematic = ((TObjString*) (*tokens)[6])->GetString();
    target = ((TObjString*) (*tokens)[5])->GetString();
    sbs_field = (((TObjString*) (*tokens)[7])->GetString()).Atoi();
    // define other things from the data map? .Atof();

    if(!(kinematic == Kin)){
      cout << "Error: The run " << run << " has a mismatch in the kinematic, investigate what is going on!" << endl;
      return;
    }

    if(!(sbs_field == SBS_field)){
      cout << "Error: The run " << run << " has a mismatch in the sbs field, investigate what is going on!" << endl;
      return;
    }

    if(!(target == targ)){
      cout << "Error: The run " << run << " has a mismatch in the target, investigate what is going on!" << endl;
      return;
    }

    if(!(pass == daPass)){
      cout << "Error: The run " << run << " has a mismatch in the pass, investigate what is going on!" << endl;
      return;
    }

    gotRun = true;

  } // end if line begins with runnum
  else{
    continue;
  }
} // end while reading lines

if((datafile.eof()) && !gotRun){
  cout << "Error:Did not find run number: " << runnum << " in the data file! Quitting, figure it out!" << endl;
  return;
}

kinematic_obj datKin(kinematic_file_name, Kin);
Ebeam = datKin.getBeamEnergy();
bbtheta = datKin.getBBAngle_Deg();
bbdist = datKin.getBBDist();
sbstheta = datKin.getSBSAngle_Deg();
sbsdist = datKin.getSBSDist();
hcaltheta = datKin.getHCalAngle_Deg();
hcaldist = datKin.getHCalDist();
Q2 = datKin.getQ2();
electron_p = datKin.getElectronP();
nucleon_p = datKin.getNucleonP();

input_file = data_object::makeInputFileName();
*/
} // end data_object

//destructor
//no dynamically allocated memory or pointers
data_object::~data_object(){}

//Implement getter functions
 int data_object::getRun(){ return run; }
 
 TString data_object::getPass(){ return pass; }
 
 TString data_object::getKinematic(){ return kinematic; }
 
 TString data_object::getTarget(){ return target; }
 
 int data_object::getSBSField(){ return sbs_field; }

 double data_object::getBeamEnergy(){ return Ebeam; }
 
 double data_object::getBBAngle_Deg(){ return bbtheta; }
 
 double data_object::getBBAngle_Rad(){ return utility::DegToRad(bbtheta); }
 
 double data_object::getBBDist(){ return bbdist; }
 
 double data_object::getSBSAngle_Deg(){ return sbstheta; }
 
 double data_object::getSBSAngle_Rad(){ return utility::DegToRad(sbstheta); }
 
 double data_object::getSBSDist(){ return sbsdist; }

 double data_object::getHCalAngle_Deg(){ return hcaltheta; }

 double data_object::getHCalAngle_Rad(){ return utility::DegToRad(hcaltheta); }

 double data_object::getHCalDist(){ return hcaldist; }

 double data_object::getQ2(){ return Q2; }

 double data_object::getElectronP(){ return electron_p; }

 double data_object::getNucleonP(){ return nucleon_p; }

 TString data_object::getInputFile(){ return input_file; }

 void data_object::printRunInfo(){
        cout << "-------------------------------------------------------------------------------------------------------------------------------------------------"        << endl
             << Form("Run number: %i,",getRun())                        << endl
             << Form("Kinematic: %s,",(getKinematic()).Data())          << endl
             << Form("Target: %s,", (getTarget()).Data())               << endl
             << Form("SBS Field: %i,",getSBSField())                    << endl
             << Form("Beam Energy: %f,",getBeamEnergy())                << endl
             << Form("BB angle in Degrees: %f,",getBBAngle_Deg())       << endl
             << Form("BB angle in Radians: %f,",getBBAngle_Rad())       << endl
             << Form("BB Distance: %f,",getBBDist())                    << endl
             << Form("SBS angle in Degrees: %f,",getSBSAngle_Deg())     << endl
             << Form("SBS angle in Radians: %f,",getSBSAngle_Rad())     << endl
             << Form("SBS Distance: %f,",getSBSDist())                  << endl
             << Form("HCal angle in Degress: %f,",getHCalAngle_Deg())   << endl
             << Form("HCal angle in Radians: %f,",getHCalAngle_Rad())   << endl
             << Form("HCal Distance: %f,",getHCalDist())                << endl
             << Form("Q2: %f,",getQ2())                                 << endl
             << Form("Electron p: %f,",getElectronP())                  << endl
             << Form("Nucleon p: %f,",getNucleonP())                    << endl
             << "-------------------------------------------------------------------------------------------------------------------------------------------------"        << endl;
 }
