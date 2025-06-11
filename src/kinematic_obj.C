//kinematic_obj.C
//Author Ezekiel Wertz
//Companion implemenation to the hearder file. The purpose is to have an object that holds all the information we care about for each kinematic

#include "../include/kinematic_obj.h"
#include "../include/utility.h"
#include <iostream>
#include "TString.h"
#include "TMath.h"

  //constructor implementation
  kinematic_obj::kinematic_obj(const char *kinematic_file_name,TString Kin){
  kinematic = Kin;
  //Now read-in the kinematic information
  ifstream kinfile(kinematic_file_name);
  //check if there is a problem opening the file
  	if(kinfile.fail()){
        cout << "Error:There was a problem with the kinematic file " << kinematic_file_name << ". Figure it out nerd!" << endl;
        return;
        }
  TString datLine;
  TString datKin;
  bool gotKin = false;
        while(datLine.ReadLine(kinfile)){
                if(datLine.BeginsWith(kinematic)){
                //We found the right kinematic. Store the info
                TObjArray *myobjs = datLine.Tokenize(",");
                //Assuming ordering is kinematic, Beam Energy, BB Angle, SBS angle, SBS dist, HCal dist, Expected Q^2, electron_p, nucleon_p
                datKin = ((TObjString*) (*myobjs)[0])->GetString();
                Ebeam = (((TObjString*) (*myobjs)[1])->GetString()).Atof();
                bbtheta = (((TObjString*) (*myobjs)[2])->GetString()).Atof();
                bbdist = (((TObjString*) (*myobjs)[3])->GetString()).Atof();
                sbstheta = (((TObjString*) (*myobjs)[4])->GetString()).Atof();
                sbsdist = (((TObjString*) (*myobjs)[5])->GetString()).Atof();
                hcaltheta = (((TObjString*) (*myobjs)[6])->GetString()).Atof();
                hcaldist = (((TObjString*) (*myobjs)[7])->GetString()).Atof();
                Q2 = (((TObjString*) (*myobjs)[8])->GetString()).Atof();
                electron_p = (((TObjString*) (*myobjs)[9])->GetString()).Atof();
                nucleon_p = (((TObjString*) (*myobjs)[10])->GetString()).Atof();
                //cout << "Kinematic " << datKin << "  Beam Energy " << Ebeam << " BB Angle  " << bbtheta << " BB Dist " << bbdist << " SBS Angle  " << sbstheta << " SBS Dist  " << sbsdist << " HCal Angle " << hcaltheta <<  " HCal Dist  " << hcaldist << " Q^2 " << Q2 << " electron_p " << electron_p << " nucelon_p " << nucelon_p  << endl;
                gotKin = true;
		}
                else{
                //Where are still searching or its a comment
                //cout << "Cond 3" << endl;
                continue;
                }//end conditional
	}//end while loop
        if ((kinfile.eof()) && !gotKin){
        //Conditional that we checked the entire kinematic file and did not find the kinematic info
        cout << "Error: Did not find kinematic: " << datKin << " in the kinematic file! Quitting, figure it out!" << endl;
        return;
        } 
  }//end constructor

  //desctructor
  //no dyanmic memory and no points
  kinematic_obj::~kinematic_obj(){}

  //getter function implementations
  TString kinematic_obj::getKinematic(){ return kinematic; }

  double kinematic_obj::getBeamEnergy(){ return Ebeam; }

  double kinematic_obj::getBBAngle_Deg(){ return bbtheta; }

  double kinematic_obj::getBBAngle_Rad(){ return utility::DegToRad(bbtheta); }

  double kinematic_obj::getBBDist(){ return bbdist; }

  double kinematic_obj::getSBSAngle_Deg(){ return sbstheta; }

  double kinematic_obj::getSBSAngle_Rad(){ return utility::DegToRad(sbstheta); }

  double kinematic_obj::getSBSDist(){ return sbsdist; }
 
  double kinematic_obj::getHCalAngle_Deg(){ return hcaltheta; }

  double kinematic_obj::getHCalAngle_Rad(){ return utility::DegToRad(hcaltheta); }

  double kinematic_obj::getHCalDist(){ return hcaldist; }

  double kinematic_obj::getQ2(){ return Q2; }
  
  double kinematic_obj::getElectronP(){ return electron_p; }
  
  double kinematic_obj::getNucleonP(){ return nucleon_p; }

  void kinematic_obj::printKinInfo(){
        cout << "-------------------------------------------------------------------------------------------------------------------------------------------------"        << endl
             << Form("Kinematic: %s,",(getKinematic()).Data())          << endl
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
