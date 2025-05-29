//exp_constants.C
//Author: Ezekiel Wertz
//The companion implementation for the header file
//Really this will just hold getter functions for any parameters that are meant to be constants but are either pass or kinematic dependent and not in config files

#include "../include/exp_constants.h"
namespace exp_constants{
  //m, height of the center of hcal above beam (m)
  const double getHCalOffset( TString pass){
  double hcal_offset;

  //First segment by pass. Later we may need to segment by kinematic
  
  	if(pass == "pass0" || pass == "pass1"){
	//If we ever want to look at this pass of data again. May need to add kinematic dependent offsets
	hcal_offset = -0.2897;
	}else if(pass == "pass2"){
	//depending on how the data looks, may need to shift this so that way all neutron peaks aline with zero
	hcal_offset = 0.0;
	}else if(pass == "mc"){
	hcal_offset = 0.0;
	}else{
 	//We should never get here, cause then we have kinematic for data we dont have
 	cout << "Error: No hcal_offset for this pass: " << pass <<" Figure it out nerd!" << endl;
	}//end conditional
  //cout << "HCal Offset: " << hcal_offset << endl;
  return hcal_offset;
  }//end function

  //maximum value of the field of the sbs magnet in Tesla. It is kinematic dependent because of magnet collision.
  const double getMaxSBSField(TString Kin, TString pass){
	double maxsbsfield = 0.0; //Tesla
	if(Kin == "SBS4" && (pass == "pass2" || pass == "mc")){
	maxsbsfield = 1.71;
	}else if(Kin == "SBS8" && pass == "pass2"){
	maxsbsfield = 1.23;
	}else if(Kin == "SBS8" && pass == "mc"){
	maxsbsfield = 1.2;
	}else if(Kin == "SBS9" && pass == "pass2"){
	maxsbsfield = 1.27;	
	}else if(Kin == "SBS9" && pass == "mc"){
	maxsbsfield = 1.24;
	}else if(( Kin == "SBS7"|| Kin == "SBS11"|| Kin == "SBS14") && (pass == "pass2" || pass == "mc")){
        maxsbsfield = 1.26;
	}else{
	//We should never get here, cause then we have kinematic for data we dont have
	cout << "Error: No sbs max field for this kinematic: " << Kin << "Figure out what you did!" << endl;
	}
  return maxsbsfield;
  }

}//end namespace
