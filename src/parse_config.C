//parse_config.C
//Author: Ezekiel Wertz
//The companion implementation for the header file
//implementation for parsing my config files. Should be able to handle all analysis types.

#include "../include/parse_config.h"
#include <iostream>

  //constructor implementation
  parse_config::parse_config(const char *setup_file_name)
  {
    ifstream setupfile(setup_file_name);
    //check if there is a problem opening the file
    if(setupfile.fail())
    {
	     cout << "ERROR: There was a problem with the setup file " << setup_file_name << ". Figure it out nerd!" << endl;
       return;
	  }
    TString myline;
    //Look for the input files. For MC we will have proton or neutron identifiers. If its data we will just have run numbers
    while(myline.ReadLine(setupfile) && !myline.BeginsWith("endrun") && !myline.BeginsWith("endfile"))
    {
	     if(!myline.BeginsWith("#")){
	     TObjArray *demObjs = myline.Tokenize(",");
	     int numObjs = demObjs->GetEntries();
       //changed to zero don't think it's an issue. Used to be 1.
		   if(numObjs>0)
       {
		       TString temp = ((TObjString*) (*demObjs)[0])->GetString();
		       bool myBool = utility::check_number(temp.Data());
		       //cout << myBool << endl;
		       //need to handle if we are dealing with MC or data
			     //We found a number
			     if(myBool)
           {
             for(int i=0;i < numObjs;i++)
             {
               //store the run numbers in the vector as ints
               int runnum_temp = (((TObjString*) (*demObjs)[i])->GetString()).Atoi();
               runnums.push_back(runnum_temp);
               //cout << runnum_temp << " ";
             }//end for loop
			     //We found a string
			     }
           else
           {
             TString key = ((TObjString*) (*demObjs)[0])->GetString();
             TString val = ((TObjString*) (*demObjs)[1])->GetString();
             if(key == "proton")
             {
               proton_root_file = val;
               //cout << "Proton File: " << proton_root_file << endl;
             }
             else if(key == "neutron")
             {
               neutron_root_file=val;
               //cout << "Neutron File: " << neutron_root_file << endl;
             }
             else if(key == "MC_file")
             {
               MC_file = val;
               //cout << "MC File: " << MC_file << endl;
             }
             else if(key == "Data_file")
             {
               Data_file = val;
             }
	     else if(key == "histfile_dir")
	     {
		histfile_dir = val;
		//cout << "Histfile_dir" << histfile_dir << endl;
	     }
             else if(key == "rootfile_dir") //keep??
             {
               rootfile_dir = val;
               //cout << "Rootfile dir" << rootfile_dir << endl;
             }
             else
             {
               //We somehow obtained a key that we were not expecting. Maybe the condition needs to be handled.
               cout << "Error: Found a key that this script can't handle. Fix that! "<< key << endl;
               return;
             }
			     }//end conditional
    		}
        else
        {
    		//We either got an empty line or 1 element.
		      cout << "Error: Line does not have the right number of elements. Look at the config file!" << endl;
    		  return;
    		}//end conditional
	     }//end conditional
    }//end of while loop
    //while loop for getting the rest of the parameters
    while(myline.ReadLine(setupfile))
    {
    	if(!myline.BeginsWith("#"))
      {
        TObjArray *daObjs = myline.Tokenize(" ");
        int nObjs = daObjs->GetEntries();
        if(nObjs >1)
        {
          TString key = ((TObjString*) (*daObjs)[0])->GetString();
          TString val = ((TObjString*) (*daObjs)[1])->GetString();
			    if(key == "exp")
          {
            Exp = val;
            //cout << "Experiment: " << Exp << endl;
          }
          else if(key == "globalcut")
          {
            globalcut = val;
			      cout << "Applying the following global cut to all data: " <<globalcut <<endl;
			    }
          else if(key == "kin")
          {
            kin = val;
            //cout << "Kinematic " << kin << endl;
          }
          else if(key == "kin_1")
          {
            kin_1 = val;
          }
          else if(key == "kin_2")
          {
            kin_2 = val;
			    }
          else if(key == "data_map_name")
          {
            data_file_name = val;
            //cout << "Data File " << data_fiel_name << endl;
          }
          else if(key == "kinematic_name")
          {
            kinematic_file_name = val;
            //cout << "Kinematic File " << kinematic_file_name << endl;
          }
          else if(key == "partial_name_p")
          {
            partial_name_p = val;
            //cout << "Partial Name P " << partial_name_p << endl;
          }
          else if(key == "partial_name_n")
          {
            partial_name_n = val;
            //cout << "Partial Name N " << partial_name_n << endl;
          }
          else if(key == "pass")
          {
            pass = val;
            //cout << "Pass " << pass << endl;
          }
          else if(key == "pass_1")
          {
            pass_1 = val;
          }
          else if(key == "pass_2")
          {
            pass_2 = val;
          }
          else if(key == "fitopt")
          {
            fitopt = val;
          }
          else if (key == "EnergyCut")
          {
            EnergyCut = val;
            //cout << "Energy Cut: " << EnergyCut << endl;
          }
          else if (key == "TrackQualityCut")
          {
            TrackQualityCut = val;
          }
          else if (key == "TrackHitsCut")
          {
            TrackHitsCut = val;
            //cout << "TrackHitsCut: " << TrackHitsCut << endl;
          }
          else if (key == "TargetVertexCut")
          {
            TargetVertexCut = val;
            //cout << "Target Vertex Cut: " << TargetVertexCut << endl;
          }
          else if (key == "W2Cut")
          {
            W2Cut = val;
            //cout << "W2 Cut: " << W2Cut << endl;
          }
          else if (key == "FidXCut")
          {
            FidXCut = val;
            //cout << "FidX Cut: " << FidXCut << endl;
          }
          else if (key == "FidYCut")
          {
            FidYCut = val;
            //cout << "FidY Cut: " << FidYCut << endl;
          }
          else if (key == "dyCut")
          {
            dyCut = val;
            //cout << "dy Cut: " << dyCut << endl;
          }
          else if (key == "eOverpCut")
          {
            eOverpCut = val;
            //cout << "eOverp Cut: " << eOverpCut << endl;
          }
          else if (key == "HCal_Energy_Cut")
          {
            HCal_Energy_Cut = val;
            //cout << "HCal_Energy_Cut: " << HCal_Energy_Cut << endl;
          }
          else if (key == "HCal_Shower_atime_Cut")
          {
            HCal_Shower_atime_Cut = val;
            //cout << "HCal_Shower_atime_Cut: " << HCal_Shower_atime_Cut << endl;
          }
          else if (key == "OpticsCut_x")
          {
            OpticsCut_x = val;
            //cout << "Optics Cut X: " << OpticsCut_x << endl;
          }
          else if (key == "OpticsCut_y")
          {
            OpticsCut_y = val;
            //cout << "Optics Cut Y: " << OpticsCut_y << endl;
          }
          else if (key == "ProtonSpotCut")
          {
            ProtonSpotCut = val;
            //cout << "ProtonSpot Cut: " << ProtonSpotCut << endl;
          }
          else if (key == "NeutronSpotCut")
          {
            NeutronSpotCut = val;
            //cout << "NeutronSpot Cut: " << NeutronSpotCut << endl;
          }
          else if (key == "isProtonCut")
          {
            isProtonCut = val;
            //cout << "isProtonCut: " << isProtonCut << endl;
          }
          else if (key == "isNeutronCut")
          {
            isNeutronCut = val;
            //cout << "isNeutronCut: " << isNeutronCut << endl;
          }
          else if(key == "HCal_Eff_map_file")
          {
            HCal_Eff_map_file = val;
            //cout << "HCal_Eff_map_file: " << HCal_Eff_map_file << endl;
          }
          else if(key == "left_right")
          {
            left_right = val;
          }
          else if(key == "SBS_field")
          {
            SBS_field = val.Atoi();
            //cout << "SBS Field " << SBS_field << endl;
          }
          else if(key == "SBS_field_1")
          {
            SBS_field_1 = val.Atoi();
          }
          else if(key == "SBS_field_2")
          {
            SBS_field_2 = val.Atoi();
          }
          else if(key == "hcalnclusmin")
          {
            hcalnclusmin = val.Atoi();
          }
          else if(key == "e_method")
          {
            e_method = val.Atoi();
            //cout << "e method " << e_method << endl;
          }
          else if(key == "slice_mode")
          {
            slice_mode = val.Atoi();
          }
          else if(key == "W2_low")
          {
            W2_low = val.Atof();
            //cout << "W2 low " << W2_low << endl;
          }
          else if(key == "W2_high")
          {
            W2_high = val.Atof();
            //cout << "W2 high " << W2_high << endl;
          }
          else if(key == "targ")
          {
            targ = val;
            //cout << "Target " << targ << endl;
          }
          else if(key == "targ_1")
          {
            targ_1 = val;
          }
          else if(key == "targ_2")
          {
            targ_2 = val;
          }
          else if(key == "spot_choice")
          {
            spot_choice = val;
          }
          else if(key == "MAXNTRACKS")
          {
            MAXNTRACKS = val.Atoi();
            //cout << "Max Number of Tracks per event " << MAXNTRACKS << endl;
          }
          else if(key == "dxO_n")
          {
            dxO_n = val.Atof();
            //cout << "x-position of neutron spot " << dxO_n << endl;
          }
          else if(key == "dyO_n")
          {
            dyO_n = val.Atof();
            //cout << "y-position of neutron spot " << dyO_n << endl;
          }
          else if(key == "dxsig_n")
          {
            dxsig_n = val.Atof();
            //cout << "x sigma of neutron spot " << dxsig_n << endl;
          }
          else if(key == "dysig_n")
          {
            dysig_n = val.Atof();
            //cout << "y sigma of neutron spot " << dysig_n << endl;
          }
          else if(key == "dxsig_fid_n")
          {
            dxsig_fid_n = val.Atof();
          }
          else if(key == "dysig_fid_n")
          {
            dysig_fid_n = val.Atof();
          }
          else if(key == "dxO_p")
          {
            dxO_p = val.Atof();
            //cout << "x-position of proton spot " << dxO_p << endl;
          }
          else if(key == "dyO_p")
          {
            dyO_p = val.Atof();
            //cout << "y-position of proton spot " << dyO_p << endl;
          }
          else if(key == "dxsig_p")
          {
            dxsig_p = val.Atof();
            //cout << "x sigma of proton spot " << dxsig_p << endl;
          }
          else if(key == "dysig_p")
          {
            dysig_p = val.Atof();
            //cout << "y sigma of proton spot " << dysig_p << endl;
          }
          else if(key == "dxsig_fid_p")
          {
            dxsig_fid_p = val.Atof();
          }
          else if(key == "dysig_fid_p")
          {
            dysig_fid_p = val.Atof();
          }
          else if(key == "dx_pn")
          {
            dx_pn = val.Atof();
            //cout << "max x difference between peaks " << dx_pn << endl;
          }
          else if(key == "useAlshield")
          {
            useAlshield = val.Atoi();
            //cout << "Use Al shield " << useAlshield << endl;
          }
          else if(key == "dx_low")
          {
            dx_low = val.Atof();
            //cout << "dx plot lower bound " << dx_low << endl;
          }
          else if(key == "dx_high")
          {
            dx_high = val.Atof();
            //cout << "dx plot higher bound " << dx_high << endl;
          }
          else if(key == "dy_low")
          {
            dy_low = val.Atof();
            //cout << "dy plot lower bound " << dy_low << endl;
          }
          else if(key == "dy_high")
          {
            dy_high = val.Atof();
            //cout << "dy plot higher bound " << dy_high << endl;
          }
          else if(key == "dxsig_n_fac")
          {
            dxsig_n_fac = val.Atof();
            //cout << "dx sigma factor for neutron " << dxsig_n_fac << endl;
          }
          else if(key == "dxsig_p_fac")
          {
            dxsig_p_fac = val.Atof();
            //cout << "dx sigma factor for proton " << dxsig_p_fac << endl;
          }
          else if(key == "dysig_n_fac")
          {
            dysig_n_fac = val.Atof();
            //cout << "dy sigma factor for neutron " << dysig_n_fac << endl;
          }
          else if(key == "dysig_p_fac")
          {
            dysig_p_fac = val.Atof();
            //cout << "dy sigma factor for proton " << dysig_p_fac << endl;
          }
          else if(key == "dysig_cut")
          {
            dysig_cut = val.Atof();
          }
          else if(key == "dysig_cut_fac")
          {
            dysig_cut_fac = val.Atof();
            //cout << "dy sigma cut factor " << dysig_cut_fac << endl;
          }
          else if(key == "spot_sig")
          {
            spot_sig = val.Atof();
            //cout << "spot sigma " << spot_sig << endl;
          }
          else if(key == "proton_thresh_fac")
          {
            proton_thresh_fac = val.Atoi();
            //cout << "Proton thresh factor: " << proton_thresh_fac << endl;
          }
          else if(key == "neutron_thresh_fac")
          {
            neutron_thresh_fac = val.Atoi();
            //cout << "Neutron thresh factor: " << neutron_thresh_fac << endl;
          }
          else if(key == "pmin")
          {
            pmin = val.Atof();
            //cout << "pmin: " << pmin << endl;
          }
          else if(key == "pmax")
          {
            pmax = val.Atof();
            //cout << "pmax: " << pmax << endl;
          }
          else if(key == "Emin")
          {
            Emin = val.Atof();
            //cout << "Emin: " << Emin << endl;
          }
          else if(key == "Emax")
          {
            Emax = val.Atof();
            //cout << "Emax: " << Emax << endl;
          }
          else if(key == "num_bin")
          {
            num_bin = val.Atoi();
            //cout << "Num bins: " << num_bin << endl;
          }
          else if(key == "coin_mean")
          {
            coin_mean = val.Atof();
            //cout << "coin mean " << coin_mean << endl;
          }
          else if(key == "coin_sigma")
          {
            coin_sigma = val.Atof();
            //cout << "coin sigma " << coin_sigma << endl;
          }
          else if(key == "coin_profile_sig")
          {
            coin_profile_sig = val.Atof();
            //cout << "coin profile sig " << coin_profile_sig << endl;
          }
          else if(key == "coin_sig_fac")
          {
            coin_sig_fac = val.Atof();
            //cout << "coin sig fac " << coin_sig_fac << endl;
          }
          else if(key == "binfac")
          {
            binfac = val.Atof();
            //cout << "binfac " << binfac << endl;
          }
          else if(key == "hbinfac")
          {
            hbinfac = val.Atof();
            //cout << "hbinfac " << hbinfac << endl;
          }
          else if(key == "W2fitmax")
          {
            W2fitmax = val.Atof();
            //cout << "W2fitmax " << W2fitmax << endl;
          }
          else if(key == "W2fitmaxwide")
          {
            W2fitmaxwide = val.Atof();
            //cout << "W2fitmaxwide " << W2fitmaxwide << endl;
          }
          else if(key == "hcalemin")
          {
            hcalemin = val.Atof();
            //cout << "hcalemin " << hcalemin << endl;
          }
          else if(key == "fidxmin")
          {
            fidx_min = val.Atof();
            //cout << "fidx_min" << fidx_min << endl;
          }
          else if(key == "fidxmax")
          {
            fidx_max = val.Atof();
            //cout << "fidx_max" << fidx_max << endl;
          }
          else if(key == "fidymin")
          {
            fidy_min = val.Atof();
            //cout << "fidy_min" << fidy_min << endl;
          }
          else if(key == "fidymax")
          {
            fidy_max = val.Atof();
            //cout << "fidy_max" << fidy_max << endl;
          }
          else if(key == "fitx_low")
          {
            fitx_low = val.Atof();
            //cout << "fitx low " << fitx_low << endl;
          }
          else if(key == "fitx_high")
          {
            fitx_high = val.Atof();
            //cout << "fitx high " << fitx_high << endl;
          }
          else if(key == "fity_low")
          {
            fity_low = val.Atof();
            //cout << "fity low " << fity_low << endl;
          }
          else if(key == "fity_high")
          {
            fity_high = val.Atof();
            //cout << "fity high " << fity_high << endl;
          }
          else if(key == "thetapq_low")
          {
            thetapq_low = val.Atof();
            //cout << "thetapq_low " << thetapq_low << endl;
          }
          else if(key == "thetapq_high")
          {
            thetapq_high = val.Atof();
            //cout << "thetapq_high " << thetapq_high << endl;
          }
          else if(key == "sf")
          {
            sf = val.Atof();
            //cout << "Scale Field " << sf << endl;
          }
          else if(key == "Ntried_override")
          {
            Ntried_override = val.Atof();
          }
          else if(key == "luminosity_override"){
            luminosity_override = val.Atof();
          }
          else if(key == "genvol_override")
          {
            genvol_override = val.Atof();
          }
          else if(key == "HCal_accep_avg_eff")
          {
            HCal_accep_avg_eff = val.Atof();
          }
          else if(key == "sync_jobs")
          {
            if(val == "true")
            {
              sync_jobs = true;
            }
            else if(val == "false")
            {
              sync_jobs = false;
            }
            else
            {
              cout << "Error: sync_jobs cannot be assigned, not boolean value!" << endl;
            }
            //cout << "Sync jobs " << sync_jobs << endl;
          }
          else if(key == "mc_override")
          {
            if(val == "true")
            {
              mc_override = true;
            }
            else if(val == "false")
            {
              mc_override = false;
            }
            else
            {
              cout << "Error: mc_override cannot be assigned, not boolean value!" << endl;
            }
          }
          else if(key == "HCal_Eff_map")
          {
            if(val == "true")
            {
              HCal_Eff_map = true;
            }
            else if(val == "false")
            {
              HCal_Eff_map = false;
            }
            else
            {
              cout << "Error: HCal_Eff_map cannot be assigned, not boolean value!" << endl;
            }
          }
          else
          {
            //We somehow obtained a key that we were not expecting. Maybe the condition needs to be handled.
            cout << "Error:Found a key that this script can't handle. Fix that! "<< key << endl;
            return;
          }
        	//remove the objects to ensure a new comes through and no runaway
        	delete daObjs;
    		}//end conditional
    		else
        {
          //We either got an empty line or 1 element.
          cout << "Line does not have the right number of elements. Look at the config file!" << endl;
          return;
    		}//end conditinal
    	}//end conditional
    }//end while

    if(runnums.empty() && (proton_root_file.Length() == 0) &&  (neutron_root_file.Length() ==0) && (rootfile_dir.Length() == 0) && (MC_file.Length() == 0) && (Data_file.Length() == 0))
    {
      // if there are data or MC files we should return true and therefore throw an error
      cout << "Error: No data files in the config file, I can't do anything if I don't know where the data lives!" << endl;
      return;
    }
  }//end constructor

  //destructor
  //no dynamic memory or pointers so no problem
  parse_config::~parse_config(){}

  //Implement get functions for each variable

  TString parse_config::getExp(){ return Exp; }

  TString parse_config::getKin(){ return kin; }

  TString parse_config::getKinFileName(){ return kinematic_file_name; }

  TString parse_config::getDataFileName(){ return data_file_name; }

  TString parse_config::getTarg(){ return targ; }

  TString parse_config::getProtonFileName(){ return proton_root_file; }

  TString parse_config::getNeutronFileName(){ return neutron_root_file; }

  TString parse_config::getMCFileName(){ return MC_file; }

  TString parse_config::getDataFile(){ return Data_file; }

  TString parse_config::getPass(){ return pass; }

  TString parse_config::getRootFileDir(){ return rootfile_dir;}

  TString parse_config::getHistFileDir(){ return histfile_dir;}

  TString parse_config::getReplayType(){ return replay_type;}

  TString parse_config::getPartialNameP(){ return partial_name_p;}

  TString parse_config::getPartialNameN(){ return partial_name_n;}

  TString parse_config::getFitOpt(){ return fitopt;}

  TString parse_config::getEnergyCut(){return EnergyCut;}

  TString parse_config::getTrackQualityCut(){return TrackQualityCut;}

  TString parse_config::getTrackHitsCut(){return TrackHitsCut;}

  TString parse_config::getTargetVertexCut(){return TargetVertexCut;}

  TString parse_config::getW2Cut(){return W2Cut;}

  TString parse_config::getFidXCut(){return FidXCut;}

  TString parse_config::getFidYCut(){return FidYCut;}

  TString parse_config::getdyCut(){return dyCut;}

  TString parse_config::geteOverpCut(){return eOverpCut;}

  TString parse_config::getHCal_Shower_atime_Cut(){return HCal_Shower_atime_Cut;}

  TString parse_config::getHCal_Energy_Cut(){return HCal_Energy_Cut;}

  TString parse_config::getOpticsCutX(){return OpticsCut_x;}

  TString parse_config::getOpticsCutY(){return OpticsCut_y;}

  TString parse_config::getProtonSpotCut(){return ProtonSpotCut;}

  TString parse_config::getNeutronSpotCut(){return NeutronSpotCut;}

  TString parse_config::getisProtonCut(){return isProtonCut;}

  TString parse_config::getisNeutronCut(){return isNeutronCut;}

  TString parse_config::get_HCalEffMapFile(){return HCal_Eff_map_file;}

  TString parse_config::get_Comp_file_1(){return Comp_file_1;}

  TString parse_config::get_Comp_file_2(){return Comp_file_2;}

  TString parse_config::get_kin_1(){return kin_1;}

  TString parse_config::get_kin_2(){return kin_2;}

  TString parse_config::get_targ_1(){return targ_1;}

  TString parse_config::get_targ_2(){return targ_2;}

  TString parse_config::get_pass_1(){return pass_1;}

  TString parse_config::get_pass_2(){return pass_2;}

  TString parse_config::get_spot_choice(){return spot_choice;}

  TString parse_config::get_left_right(){return left_right;}

  int parse_config::getSBSField(){ return SBS_field; }

  int parse_config::getAlshield(){ return useAlshield; }

  int parse_config::getMAXNTRACKS(){ return MAXNTRACKS; }

  int parse_config::get_emethod(){ return e_method; }

  int parse_config::get_HCalNclusMin(){return hcalnclusmin; }

  int parse_config::getSBS_field_1(){return SBS_field_1;}

  int parse_config::getSBS_field_2(){return SBS_field_2;}

  int parse_config::get_slice_mode(){return slice_mode;}

  double parse_config::get_dxOn(){ return dxO_n; }

  double parse_config::get_dyOn(){ return dyO_n; }

  double parse_config::get_dxsign(){ return dxsig_n; }

  double parse_config::get_dysign(){ return dysig_n; }

  double parse_config::get_dxsigfidn(){return dxsig_fid_n;}

  double parse_config::get_dysigfidn(){return dysig_fid_n;}

  double parse_config::get_dxOp(){ return dxO_p; }

  double parse_config::get_dyOp(){ return dyO_p; }

  double parse_config::get_dxsigp(){ return dxsig_p; }

  double parse_config::get_dysigp(){ return dysig_p; }

  double parse_config::get_dxsigfidp(){return dxsig_fid_p;}

  double parse_config::get_dysigfidp(){return dysig_fid_p;}

  double parse_config::get_dxpn(){ return dx_pn; }

  double parse_config::getW2Low(){ return W2_low; }

  double parse_config::getW2High(){ return W2_high; }

  double parse_config::get_dxLow(){ return dx_low; }

  double parse_config::get_dxHigh(){ return dx_high; }

  double parse_config::get_dyLow(){ return dy_low; }

  double parse_config::get_dyHigh(){ return dy_high; }

  double parse_config::get_dxSignFac(){ return dxsig_n_fac; }

  double parse_config::get_dxSigpFac(){ return dxsig_p_fac; }

  double parse_config::get_dySignFac(){ return dysig_n_fac; }

  double parse_config::get_dySigpFac(){ return dysig_p_fac; }

  double parse_config::get_dySigCut(){ return dysig_cut; }

  double parse_config::get_dySigCutFac(){ return dysig_cut_fac; }

  double parse_config::getCoinMean(){ return coin_mean; }

  double parse_config::getCoinSig(){ return coin_sigma; }

  double parse_config::getCoinProfSig(){ return coin_profile_sig; }

  double parse_config::getCoinSigFac(){ return coin_sig_fac; }

  double parse_config::getHCaleMin(){ return hcalemin; }

  double parse_config::getProtonThreshFac(){ return proton_thresh_fac; }

  double parse_config::getNeutronThreshFac(){ return neutron_thresh_fac; }

  double parse_config::getNumBin(){ return num_bin; }

  double parse_config::getPmax(){ return pmax; }

  double parse_config::getPmin(){ return pmin; }

  double parse_config::getEmin(){ return Emin; }

  double parse_config::getEmax(){ return Emax; }

  double parse_config::getThetapqLow(){ return thetapq_low; }

  double parse_config::getThetapqHigh(){ return thetapq_high; }

  double parse_config::getW2FitMax(){ return W2fitmax; }

  double parse_config::getW2FitMaxWide(){ return W2fitmaxwide; }

  double parse_config::getBinFac(){ return binfac; }

  double parse_config::getHBinFac(){ return hbinfac; }

  double parse_config::get_sf(){ return sf; }

  double parse_config::getNTriedOverride(){return Ntried_override;}

  double parse_config::getLumiOverride(){return luminosity_override;}

  double parse_config::getVolOverride(){return genvol_override;}

  double parse_config::get_HCalAccepAvgEff(){return HCal_accep_avg_eff;}

  double parse_config::get_fidxmin(){return fidx_min;}

  double parse_config::get_fidxmax(){return fidx_max;}

  double parse_config::get_fidymin(){return fidy_min;}

  double parse_config::get_fidymax(){return fidy_max;}

  double parse_config::get_fitxlow(){return fitx_low;}

  double parse_config::get_fitxhigh(){return fitx_high;}

  double parse_config::get_fitylow(){return fity_low;}

  double parse_config::get_fityhigh(){return fity_high;}

  double parse_config::get_spotsig(){return spot_sig;}

  bool parse_config::get_MCOverride(){return mc_override;}

  bool parse_config::get_syncJobs(){ return sync_jobs; }

  bool parse_config::get_HCalEffMap(){return HCal_Eff_map;}

  TCut parse_config::getGlobalCut(){ return globalcut; }

  vector<int> parse_config::getRunNums(){ return runnums; }

  void parse_config::printRunNums()
  {
  	for(int j=0; j< runnums.size(); j++)
    {
      if(j==runnums.size()-1)
      {
        cout << runnums[j] << endl;
      }
      else
      {
        cout << runnums[j] << ", ";
      }
  	}
  }
