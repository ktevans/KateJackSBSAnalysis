#ifndef DBparse_H
#define DBparse_H

# Based off of a script from Sean Jeffas. A lot of the information references here is not useful, but it helps with He3 polarization.

namespace DBparse {

  TString DB_dir = "/w/halla-scshelf2102/sbs/ktevans/KateJackSBSAnalysis/config/";
  TString DB_corr_dir = "/w/halla-scshelf2102/sbs/ktevans/KateJackSBSAnalysis/config/corrections/";
  
  std::map<TString, TString> DBFileMap {
    {"He3 Polarization", "He3_pol.csv"},
    {"Beam Polarization", "Beam_pol.csv"},
    {"Helicity Quality", "Helicity_quality.csv"},
    {"Moller Quality", "Moller_quality.csv"},
    {"Field Measurement", "Field_Meas.csv"},
    {"Asymmetry Correction", "corr"}
  };
  
  struct DBrequest{
    TString var_names;   // Variable name
    TString info;        // More description about this variable
    bool    mandatory;   // Is variable mandatory?
  };

  struct DBInfo{
    TString                  cfg;
    double                   W2min;
    double                   W2max;
    double                   dymax;
    vector<DBrequest>        var_req;
    map<TDatime,double>      He3Pol;
    vector<vector<TDatime*>> BeamPolTime;
    vector<vector<double>>   BeamPolValue;
    map<int,double>          GoodHel;
    map<int,double>          GoodMoller;
    double                   AccidentalAsymmetry = 0;
    double                   AccidentalAsymmetryErr = 0;
    double                   AccidentalFraction = 0;
    double                   AccidentalFractionErr = 0;
    double                   PionAsymmetry = 0;
    double                   PionAsymmetryErr = 0;
    double                   PionFraction = 0;
    double                   PionFractionErr = 0;
    double                   InelasticAsymmetry = 0;
    double                   InelasticAsymmetryErr = 0;
    double                   InelasticFraction = 0;
    double                   InelasticFractionErr = 0;
    double                   NitrogenFraction = 0;
    double                   NitrogenFractionErr = 0;
  };


  void DB_load(DBInfo &request);
  void DB_SetCorrections(DBInfo &DBInfo);
 
}

#endif
