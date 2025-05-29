//data_object.C
//Author: Ezekiel Wertz
//Companion implementation. Class to represent a data object file for GMn/nTPE analysis. Anything we need to know about data files.

#include "../include/data_object.h"
#include <iostream>
#include "TMath.h"


//private helper function. Requires that private class variables are initialized first
//need to test this function carefully. Should be correct for data file directory structure
TString data_object::makeInputFileName(){
 //All of this was to have a modular input directory. So let's make it
 string input_directory_1 = "/work/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take3";
 string input_directory_2 = "/work/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take4";
 TString inputfile;
 const char *input_directory_1_char = input_directory_1.c_str();
 const char *input_directory_2_char = input_directory_2.c_str();
 const char *pass_char = pass.Data();
 const char *kin_char = kinematic.Data();
 const char *tar_char = target.Data();
 	if((pass == "pass0") || (pass == "pass1")){
	//pass0 and pass1 data is no longer on work. You have to bring it back from cache yourself
	cout << "Pass0 and Pass1 no longer on work. Bring it back from tape!" << endl;
	//This file structure is still correct so I'm not gonna get rid of it
 	inputfile = Form("/work/halla/sbs/sbs-gmn/%s/%s/%s/rootfiles/e1209019_fullreplay_%i_*.root",pass_char,kin_char,tar_char,run);
	}else if(pass == "pass2"){
	//need to implement special treatment for all kinematics. Due to different file structure.
		//Only implemented for files in the LH2 and LD2 directories
		if(kinematic == "SBS8" && ((target == "LH2")||(target == "LD2"))){
			//Then its segmented by magnet field setting
			//Then handle SBS 70%
			
			//SBS 70%
			if(sbs_field == 70 ){
			
				//now we need to handle the different parts in SBS 70%
				//The best way to do this is probably by run number
				if((run >= 13444) && (run <= 13455)){
				inputfile = Form("%s/%s/%s/rootfiles/SBS%ipercent_part1/e1209019_fullreplay_%i_*.root",input_directory_1_char,kin_char,tar_char,sbs_field,run);
				}else if((run >= 13482) && (run <= 13505)){
				inputfile = Form("%s/%s/%s/rootfiles/SBS%ipercent_part2/e1209019_fullreplay_%i_*.root",input_directory_1_char,kin_char,tar_char,sbs_field,run);
				}else if((run >= 13558) && (run <= 13578)){
				inputfile = Form("%s/%s/%s/rootfiles/SBS%ipercent_part3/e1209019_fullreplay_%i_*.root",input_directory_1_char,kin_char,tar_char,sbs_field,run);
				}else if((run >= 13587) && (run <= 13620)){
				inputfile = Form("%s/%s/%s/rootfiles/SBS%ipercent_part4/e1209019_fullreplay_%i_*.root",input_directory_1_char,kin_char,tar_char,sbs_field,run);
				}else{
				//We probably got here only if we found a run number that does not meet the other conditionals. So let's give it an error message and provide a default input file similar to the other part of the conditonal
				cout << "Error: Got a run number that is not in any range for SBS8 70% " << kinematic << " " << sbs_field << " " << run << endl;
				inputfile = Form("%s/%s/%s/rootfiles/SBS%ipercent/e1209019_fullreplay_%i_*.root",input_directory_1_char,kin_char,tar_char,sbs_field,run);
					
				}//end conditional
			//Handles SBS 0%, 50%, 100%
			}else if((sbs_field == 0)||(sbs_field == 50)||(sbs_field == 100)){
			inputfile = Form("%s/%s/%s/rootfiles/SBS%ipercent/e1209019_fullreplay_%i_*.root",input_directory_1_char,kin_char,tar_char,sbs_field,run);			
			}else{
			//make a catch all error if for some reason we get a magnetic field that doesnt make sense
			cout << "Error: Trouble with finding the data file! Kinematic: " << kinematic << " Target: " << target << " SBS Field: " << sbs_field << endl;
			}//end conditional
		//need to handle SBS8 dummy and optics properly
		}else if(kinematic == "SBS8" && ((target == "Dummy")||(target == "Optics"))){
		inputfile = Form("%s/%s/%s/rootfiles/e1209019_fullreplay_%i_*.root",input_directory_1_char,kin_char,tar_char,run);
		//handle everything for SBS7 and SBS9
		}else if(kinematic == "SBS7" || kinematic == "SBS9"){
		inputfile = Form("%s/%s/%s/rootfiles/e1209019_fullreplay_%i_*.root",input_directory_1_char,kin_char,tar_char,run);
		//Need to handle LD2 for SBS11
		}else if(kinematic == "SBS11" && (target == "LD2")){
		//The best way to do this is probably by run number
			if((run >= 12314) && (run <= 12425)){
                	inputfile = Form("%s/%s/%s/part1/rootfiles/e1209019_fullreplay_%i_*.root",input_directory_2_char,kin_char,tar_char,run);
               		}else if((run >= 12473) && (run <= 12830)){
               		//really silly file structure here but we must handle it
				if((run >= 12473) && (run <= 12499)){
				inputfile = Form("%s/%s/%s/part2/run124XX/rootfiles/e1209019_fullreplay_%i_*.root",input_directory_2_char,kin_char,tar_char,run);
				}else if((run >= 12500) && (run <= 12599)){
				inputfile = Form("%s/%s/%s/part2/run125XX/rootfiles/e1209019_fullreplay_%i_*.root",input_directory_2_char,kin_char,tar_char,run);
				}else if((run >= 12614) && (run <= 12699)){
                                inputfile = Form("%s/%s/%s/part2/run126XX/rootfiles/e1209019_fullreplay_%i_*.root",input_directory_2_char,kin_char,tar_char,run);
                                }else if((run >= 12700) && (run <= 12799)){
                                inputfile = Form("%s/%s/%s/part2/run127XX/rootfiles/e1209019_fullreplay_%i_*.root",input_directory_2_char,kin_char,tar_char,run);
                                }else if((run >= 12801) && (run <= 12830)){
                                inputfile = Form("%s/%s/%s/part2/run128XX/rootfiles/e1209019_fullreplay_%i_*.root",input_directory_2_char,kin_char,tar_char,run);
                                }else{
				cout << "Error: We got a run number: " << run << " that we cannot handle for SBS11 LD2!" << endl;
				}
				
			}else if((run >= 12894) && (run <= 13063)){
                        inputfile = Form("%s/%s/%s/part3/rootfiles/e1209019_fullreplay_%i_*.root",input_directory_2_char,kin_char,tar_char,run);
               		}else{
			cout << "Error: Got a run number that is not in any range for SBS11 LD2 " << kinematic << " " << sbs_field << " " << run << endl;
			inputfile = Form("%s/%s/%s/rootfiles/e1209019_fullreplay_%i_*.root",input_directory_2_char,kin_char,tar_char,run);
			}//end conditional
		//need to handle all other targets SBS11
		}else if(kinematic == "SBS11" && (target == "LH2" || target == "Optics" || target == "Dummy")){	
		inputfile = Form("%s/%s/%s/rootfiles/e1209019_fullreplay_%i_*.root",input_directory_2_char,kin_char,tar_char,run);
		//handle SBS14 and SBS4
		}else if(kinematic == "SBS14" || kinematic == "SBS4"){
		inputfile = Form("%s/%s/%s/rootfiles/e1209019_fullreplay_%i_*.root",input_directory_2_char,kin_char,tar_char,run);
		}else{
		//We should never get here unless there is an error, all above conditionals handle file structure already since it is not great.
			cout << "Error: Got to file making function else statement, something is most likely wrong." << endl; 
			if(kinematic == "SBS7" || kinematic == "SBS8" || kinematic == "SBS9"){
			inputfile = Form("%s/%s/%s/rootfiles/e1209019_fullreplay_%i_*.root",input_directory_1_char,kin_char,tar_char,run);
			}else if(kinematic == "SBS11" || kinematic == "SBS14" || kinematic == "SBS4"){
			inputfile = Form("%s/%s/%s/rootfiles/e1209019_fullreplay_%i_*.root",input_directory_2_char,kin_char,tar_char,run);
			}else {
			//We should never get here
			cout << "Error: Something is serioulsy wrong! Cannot find your data file." << endl;
			}
		}//end kinematic conditionals
	}else{
	//make an error. Some how we got not pass 0,1, or 2
	cout << "Error: Pass variable was given that is not pass 0,1, or 2 " << pass << " !" << endl;
	}
 //cout << "File Location " << inputfile << endl;
 return inputfile;
}//end input file maker

//constructor implementation
data_object::data_object(int runnum,const char *data_file_name,const char *kinematic_file_name,TString Kin, TString SBS_field, TString targ, TString daPass){
 //This part of the constructor reads in information about the data file itself
 ifstream datafile(data_file_name);
 //check if there is a problem opening the file
 	if(datafile.fail()){
 	cout << "Error:There was a problem with the data file " << data_file_name << ". Figure it out nerd!" << endl;
 	return;
 	}
 TString currentLine;
 bool gotRun = false;
 TString runnum_string = utility::intToTString(runnum);
 	while(currentLine.ReadLine(datafile)){
        	if(currentLine.BeginsWith(runnum_string)){
		//We found the right run number. Store the info
		TObjArray *tokens = currentLine.Tokenize(" ");
        	//Assuming ordering is runnum, pass, SBSKin, Target, sbsfield
		run =  (((TObjString*) (*tokens)[0])->GetString()).Atoi();
         	pass = ((TObjString*) (*tokens)[1])->GetString();
         	kinematic =  ((TObjString*) (*tokens)[2])->GetString();
         	target =  ((TObjString*) (*tokens)[3])->GetString();
         	sbs_field = (((TObjString*) (*tokens)[4])->GetString()).Atoi();
         	//cout << "Run " << run  << " Pass " << pass << " Kin " << kinematic << " Target " << target << " SBS Field  " << sbs_field << endl;
			if(!(kinematic == Kin)){
        		cout << "Error: The run " << run << " has a mismatch in the kinematic, investigate what is going on!" << endl;
        		return;
        		}
        		if(!(sbs_field == SBS_field)){
        		cout << "Error: The run "<< run << " has a mismatch in the sbs field, investigate what is going on!" << endl;
        		return;
        		}
        		if(!(target == targ )){
        		cout << "Error: The run " << run << " has a mismatch in the target, investigate what is going on!" << endl;
        		return;
        		}
			if(!(pass == daPass )){
			cout << "Error: The run " << run << " has a mismatch in the pass, investigate what is going on!" << endl;
			return;
			}
        	gotRun = true;
         	//cout << gotRun << endl;
         	} else{
         	//We are still searching but it's a comment denoted by #
         	//cout << "Cond 2" << endl; 
         	continue;
         	}//end conditional
	}//end while loop
 if ((datafile.eof()) && !gotRun){
 //Conditional that we checked the entire data file and did not find the runnum
 cout << "Error:Did not find run number: " << runnum << " in the data file! Quitting, figure it out!" << endl;
 return;
 }
 //Now read-in the kinematic information and store the info
 kinematic_obj datKin(kinematic_file_name, Kin);
 Ebeam = datKin.getBeamEnergy();
 bbtheta = datKin.getBBAngle_Deg();
 bbdist = datKin.getBBDist();
 sbstheta = datKin.getSBSAngle_Deg();
 sbsdist = datKin.getSBSDist();
 hcaltheta = datKin.getHCalAngle_Deg();
 hcaldist = datKin.getHCalDist();
 Q2=datKin.getQ2();
 electron_p = datKin.getElectronP();
 nucleon_p = datKin.getNucleonP();

 input_file = data_object::makeInputFileName();

}//end constructor

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
