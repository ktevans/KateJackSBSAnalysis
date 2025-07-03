#ifndef DBparse_H
#define DBparse_H

// Based off of a script from Sean Jeffas. A lot of the information references here is not useful, but it helps with He3 polarization.

namespace DBparse {

  TString DB_dir = "/w/halla-scshelf2102/sbs/ktevans/KateJackSBSAnalysis/config/";
  //TString DB_corr_dir = "/w/halla-scshelf2102/sbs/ktevans/KateJackSBSAnalysis/config/corrections/";
  
  std::map<TString, TString> DBFileMap {
    {"He3 Polarization", "He3_pol.csv"}
  };

//{"Asymmetry Correction", "corr"}, {"Field Measurement", "Field_Meas.csv"}, {"Helicity Quality", "Helicity_quality.csv"}, {"Moller Quality", "Moller_quality.csv"}, {"Beam Polarization", "Beam_pol.csv"}
  
  struct DBrequest{
    TString var_names;   // Variable name
    TString info;        // More description about this variable
    bool    mandatory;   // Is variable mandatory?
  };

  struct DBInfo{
    TString                  cfg;
    vector<DBrequest>        var_req;
    map<TDatime,double>      He3Pol;
  };


  void DB_load(DBInfo &request);
  void DB_SetCorrections(DBInfo &DBInfo);
 
}

#endif
