//Original Author Ezekiel Wertz
//Update Author Kate Evans
//05/29/2025
//Purpose: Parsing Script for GEn-II data to produce output histograms for later analysis

//The exact ordering of this matters. ROOT for some reason cannot handle calls for files that have already been included.
#include "TF1.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCut.h"
#include "TEventList.h"
#include "TStopwatch.h"
#include "TTreeFormula.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "../src/utility.C"
#include "../src/exp_constants.C"
#include "../src/kinematic_obj.C"
#include "../src/data_object.C"
#include "../src/cuts.C"
#include "../src/physics.C"
#include "../src/parse_config.C"
#include "../include/DBparse.h"

DBparse::DBInfo DBInfo;

// Load database files
void getDB(TString cfg){
  
  cout<<"Attempting to load DB File"<<endl;
  cout<<"---------------------------------------------------------------"<<endl;

   vector<DBparse::DBrequest> request = {
     {"He3 Polarization","He3 target polarization", 1},
  };

  DBInfo.cfg = cfg;
  DBInfo.var_req = request;

  DB_load(DBInfo);

  cout<<"---------------------------------------------------------------"<<endl;

}

//Main
void data_parse(const char *setup_file_name){

  //Define a clock to check macro processing time
  TStopwatch *watch = new TStopwatch();
  watch->Start( kTRUE );

  //parse object to get in the information that The One Config file has and is manipulated
  parse_config mainConfig(setup_file_name);
  //mainConfig.printDataYields();

  //store all the parameters from the mainConfig file into local variables. So we don't have to keep recalling them
  vector<int> runNums = mainConfig.getRunNums();
  TCut globalcut = mainConfig.getGlobalCut();
  TString exp = mainConfig.getExp();
  TString kin = mainConfig.getKin();
  TString kinematic_file = mainConfig.getKinFileName();
  TString data_map = mainConfig.getDataFileName();
  TString pass = mainConfig.getPass();
  int sbs_field = mainConfig.getSBSField();
  int maxtracks = mainConfig.getMAXNTRACKS();
  TString target = mainConfig.getTarg();
  double W2_low = mainConfig.getW2Low();
  double W2_high = mainConfig.getW2High();
  double dxO_n = mainConfig.get_dxOn();
  double dyO_n = mainConfig.get_dyOn();
  double dxsig_n = mainConfig.get_dxsign();
  double dysig_n = mainConfig.get_dysign();
  double dxO_p = mainConfig.get_dxOp();
  double dyO_p = mainConfig.get_dyOp();
  double dxsig_p = mainConfig.get_dxsigp();
  double dysig_p = mainConfig.get_dysigp();
  double dx_low = mainConfig.get_dxLow();
  double dx_high = mainConfig.get_dxHigh();
  double dy_low = mainConfig.get_dyLow();
  double dy_high = mainConfig.get_dyHigh();
  int useAlshield = mainConfig.getAlshield();
  double dxsig_n_fac = mainConfig.get_dxSignFac();
  double dxsig_p_fac = mainConfig.get_dxSigpFac();
  double dysig_n_fac = mainConfig.get_dySignFac();
  double dysig_p_fac = mainConfig.get_dySigpFac();
  double W2fitmax = mainConfig.getW2FitMax();
  double binfac = mainConfig.getBinFac();
  double hbinfac = mainConfig.getHBinFac();
  double hcal_offset = exp_constants::getHCalOffset(pass);
  int e_method = mainConfig.get_emethod();
  double coin_mean = mainConfig.getCoinMean();
  double coin_sig_fac = mainConfig.getCoinSigFac();
  double coin_profile_sig = mainConfig.getCoinProfSig();
  double coin_sig =  mainConfig.getCoinSig();
  double hcalemin = mainConfig.getHCaleMin();
  double dysig_cut = mainConfig.get_dySigCut();
  double dysig_cut_fac = mainConfig.get_dySigCutFac();
  int hcalnclusmin = mainConfig.get_HCalNclusMin();


  //store all important kinematic info in local variables
  kinematic_obj myKin(kinematic_file,kin);
  double EBeam = myKin.getBeamEnergy();
  double hcaldist = myKin.getHCalDist();
  double sbsdist = myKin.getSBSDist();
  double bbtheta = myKin.getBBAngle_Rad();
  double hcaltheta = myKin.getHCalAngle_Rad();
  double p_nuc_centr = myKin.getNucleonP();

  //setup hcal physical bounds that match database for each pass
  vector<double> hcalpos = cuts::hcal_Position_data(pass);

  //setup hcal active area with bounds that match database depending on pass
  vector<double> hcalaa = cuts::hcal_ActiveArea_data(1,1,pass);

  //setup fiducial region based on dx and dy spot information
  vector<double> hcalfid = cuts::hcalfid(dxsig_p,dxsig_n,dysig_p,hcalaa,dxsig_p_fac,dysig_p_fac);

  //print out the information for the fid region
  for(int k=0; k<hcalfid.size();k++){
  cout <<"HCal Fid "<< k << " :" << hcalfid[k] << endl;
  }

  int num_runs = runNums.size();
  //double hcalfit_low = exp_constants::hcalposXi_mc; //lower fit/bin limit for hcal dx plots.
  //double hcalfit_high = exp_constants::hcalposXf_mc; //higher fit/bin limit for hcal dx plots.
  //double hcal_fitrange = exp_constants::hcal_vrange; //Full range of hcal dx plots

  //store all the data information we care about in a vector of data_objects for further use
  vector<data_object> myData;
  for(int i=0; i < num_runs; i++ ){
  //Run Num, date map name, kinematic map name, Kinematic, SBS Field, Target, Pass
  data_object myObj(runNums[i],data_map,kinematic_file,kin,sbs_field,target,pass);
  myData.push_back(myObj);
  myObj.printRunInfo();
  }

  //setup output file
  TString outfile = utility::makeOutputFileNameParse(exp,pass,kin,sbs_field,target);
  TFile *fout = new TFile(outfile,"RECREATE");

  double dx_pn = mainConfig.get_dxpn();

  //allocate memory at each run
  TChain *C = nullptr;

  //create output tree
  TTree *Parse = new TTree("Parse","Analysis Tree");

  //uncut output tree variables

  //calculated and CODA variables
  double dx_out, dy_out, xexp_out, yexp_out, W2_out, nu_out, tau_out, epsilon_out, pcorr_out, mott_out;
  Parse->Branch("dx", &dx_out, "dx/D");
  Parse->Branch("dy", &dy_out, "dy/D");
  Parse->Branch("sbs.hcal.x_exp", &xexp_out, "sbs.hcal.x_exp/D");
  Parse->Branch("sbs.hcal.y_exp", &yexp_out, "sbs.hcal.y_exp/D");
  Parse->Branch("e.kine.W2", &W2_out, "e.kine.W2/D");
  Parse->Branch("e.kine.nu", &nu_out, "e.kine.nu/D");
  Parse->Branch("e.kine.tau", &tau_out, "e.kine.tau/D");
  Parse->Branch("e.kine.epsilon", &epsilon_out, "e.kine.epsilon/D");
  Parse->Branch("pcorr", &pcorr_out, "pcorr/D");
  Parse->Branch("e.kine.mott", &mott_out, "e.kine.mott/D");

  double Ebeam_corr_out, E_eprime_out, etheta_out, p_N_out, Q2_out;
  Parse->Branch("Ebeam_corr", &Ebeam_corr_out, "Ebeam_corr/D");
  Parse->Branch("E_eprime", &E_eprime_out, "E_eprime/D");
  Parse->Branch("etheta", &etheta_out, "etheta/D");
  Parse->Branch("p_N", &p_N_out, "p_N/D");
  Parse->Branch("e.kine.Q2", &Q2_out, "e.kine.Q2/D");

  double He3Pol_out, g_trigbits_out, g_evtime_out;
  int helicity_out;
  Parse->Branch("he3pol", &He3Pol_out, "he3pol/D");
  Parse->Branch("g.trigbits", &g_trigbits_out, "g.trigbits/D");
  Parse->Branch("g.evtime", &g_evtime_out, "g.evtime/D");
  Parse->Branch("helicity", &helicity_out, "helicity/I");

  TDatime datetime_out;
  Parse->Branch("datetime", "TDatime", &datetime_out);


  //hcal variables
  double sbs_hcal_x_out, sbs_hcal_y_out, sbs_hcal_e_out, sbs_hcal_atimeblk_out, sbs_hcal_clus_blk_row_out, sbs_hcal_clus_blk_col_out, sbs_hcal_nblk_out;//,ehcal_tree_out;
  Parse->Branch("sbs.hcal.x", &sbs_hcal_x_out, "sbs.hcal.x/D");
  Parse->Branch("sbs.hcal.y", &sbs_hcal_y_out, "sbs.hcal.y/D");
  Parse->Branch("sbs.hcal.e", &sbs_hcal_e_out, "sbs.hcal.e/D");
  //Parse->Branch("ehcal_tree", &ehcal_tree_out, "ehcal_tree/D");
  Parse->Branch("sbs.hcal.atimeblk", &sbs_hcal_atimeblk_out, "sbs.hcal.atimeblk/D");
  Parse->Branch("sbs.hcal.nblk", &sbs_hcal_nblk_out, "sbs.hcal.nblk/D");
  Parse->Branch("sbs.hcal.clus_blk.row",&sbs_hcal_clus_blk_row_out, "sbs.hcal.clus_blk.row/D");
  Parse->Branch("sbs.hcal.clus_blk.col", &sbs_hcal_clus_blk_col_out, "sbs.hcal.clus_blk.col/D");

  int Ndata_sbs_hcal_clus_id_out, sbs_hcal_nclus_out, sbs_hcal_clus_blk_id_out;
  Parse->Branch("Ndata.sbs.hcal.clus.id", &Ndata_sbs_hcal_clus_id_out, "Ndata.sbs.hcal.clus.id/I");
  Parse->Branch("sbs.hcal.clus_blk.id", &sbs_hcal_clus_blk_id_out, "sbs.hcal.clus_blk.id/I");
  Parse->Branch("sbs.hcal.nclus", &sbs_hcal_nclus_out,"sbs.hcal.nclus/I");


  //bbcal variables
  double BBtot_e_out, bb_sh_e_out, bb_ps_e_out, bb_sh_atimeblk_out, bb_ps_atimeblk_out;
  Parse->Branch("bb.tot.e", &BBtot_e_out, "bb.tot.e/D");
  Parse->Branch("bb.sh.e", &bb_sh_e_out, "bb.sh.e/D");
  Parse->Branch("bb.ps.e", &bb_ps_e_out, "bb.ps.e/D");
  Parse->Branch("bb.sh.atimeblk", &bb_sh_atimeblk_out, "bb.sh.atimeblk/D");
  Parse->Branch("bb.ps.atimeblk", &bb_ps_atimeblk_out, "bb.ps.atimeblk/D");

  int bb_sh_nclus_out, bb_sh_nblk_out, bb_ps_nclus_out, bb_ps_nblk_out;
  Parse->Branch("bb.sh.nclus", &bb_sh_nclus_out, "bb.sh.nclus/I");
  Parse->Branch("bb.sh.nblk", &bb_sh_nblk_out, "bb.sh.nblk/I");
  Parse->Branch("bb.ps.nclus", &bb_ps_nclus_out, "bb.ps.nclus/I");
  Parse->Branch("bb.ps.nblk", &bb_ps_nblk_out, "bb.ps.nblk/I");

  //grinch variables
  double bb_grinch_tdc_clus_trackindex_out, bb_grinch_tdc_clus_size_out, bb_grinch_tdc_hit_time_out[maxtracks], bb_grinch_tdc_hit_pmtnum_out[maxtracks], bb_grinch_tdc_hit_amp_out[maxtracks];
  Parse->Branch("bb.grinch_tdc.clus.trackindex",&bb_grinch_tdc_clus_trackindex_out,"bb.grinch_tdc.clus.trackindex/D");
  Parse->Branch("bb.grinch_tdc.clus.size",&bb_grinch_tdc_clus_size_out,"bb.grinch_tdc.clus.size/D");
  Parse->Branch("bb.grinch_tdc.hit.time",&bb_grinch_tdc_hit_time_out,"bb.grinch_tdc.hit.time/D");
  Parse->Branch("bb.grinch_tdc.hit.pmtnum",&bb_grinch_tdc_hit_pmtnum_out,"bb.grinch_tdc.hit.pmtnum/D");
  Parse->Branch("bb.grinch_tdc.hit.amp",&bb_grinch_tdc_hit_amp_out,"bb.grinch_tdc.hit.amp/D");
  int Ndata_bb_grinch_tdc_hit_time_out;
  Parse->Branch("Ndata.bb.grinch_tdc.hit.time",&Ndata_bb_grinch_tdc_hit_time_out,"Ndata.bb.grinch_tdc.hit.time/I");


  //gem variables
  double bb_gem_track_nhits_out, bb_gem_track_ngoodhits_out, bb_gem_track_chi2ndf_out;
  Parse->Branch("bb.gem.track.nhits", &bb_gem_track_nhits_out, "bb.gem.track.nhits/D");
  Parse->Branch("bb.gem.track.ngoodhits", &bb_gem_track_ngoodhits_out, "bb.gem.track.ngoodhits/D");
  Parse->Branch("bb.gem.track.chi2ndf", &bb_gem_track_chi2ndf_out, "bb.gem.track.chi2ndf/D");


  //tracking variables
  double bb_tr_x_out, bb_tr_y_out, bb_tr_p_out, bb_tr_vz_out, bb_tr_th_out, bb_tr_ph_out, bb_tr_r_x_out, bb_tr_r_y_out, bb_tr_r_th_out, bb_tr_r_ph_out;
  Parse->Branch("bb.tr.x", &bb_tr_x_out, "bb.tr.x/D");
  Parse->Branch("bb.tr.y", &bb_tr_y_out, "bb.tr.y/D");
  Parse->Branch("bb.tr.p", &bb_tr_p_out, "bb.tr.p/D");
  Parse->Branch("bb.tr.vz", &bb_tr_vz_out, "bb.tr.vz/D");
  Parse->Branch("bb.tr.th", &bb_tr_th_out, "bb.tr.th/D");
  Parse->Branch("bb.tr.ph", &bb_tr_ph_out, "bb.tr.ph/D");
  Parse->Branch("bb.tr.r_x", &bb_tr_r_x_out, "bb.tr.r_x/D");
  Parse->Branch("bb.tr.r_y", &bb_tr_r_y_out, "bb.tr.r_y/D");
  Parse->Branch("bb.tr.r_th", &bb_tr_r_th_out, "bb.tr.r_th/D");
  Parse->Branch("bb.tr.r_ph", &bb_tr_r_ph_out, "bb.tr.r_ph/D");

  int bb_tr_n_out;
  Parse->Branch("bb.tr.n", &bb_tr_n_out, "bb.tr.n/I");

  double BB_E_over_p_out;
  Parse->Branch("bb.e_over_p", &BB_E_over_p_out, "bb.e_over_p/D");

  //kinematic variables
  int run_out, mag_out;
  Parse->Branch("runnum", &run_out, "runnum/I");
  Parse->Branch("sbs.field", &mag_out, "sbs.field/I");

  double proton_deflection_out, p_central_out;
  Parse->Branch("proton_deflection", &proton_deflection_out, "proton_deflection/D");
  Parse->Branch("p_central", &p_central_out, "p_centrial/D");

  //cut variables
  int passGlobal_out, HCalON_out;
  Parse->Branch("passGlobal", &passGlobal_out,"passGlobal/I");
  Parse->Branch("HCalON", &HCalON_out,"HCalON/I");

  double hcal_sh_atime_diff_out;
  Parse->Branch("adc.coin", &hcal_sh_atime_diff_out, "adc.coin/D");


  //double coin_mean_out;
  //double coin_sigma_out;
  //double dyO_p_out;
  //double dyO_n_out;
  //double dysig_p_out;
  //double dysig_n_out;
  //double dysig_n_fac_out;
  //double dysig_p_fac_out;
  //double dxO_p_out;
  //double dxO_n_out;
  //double dxsig_p_out;
  //double dxsig_n_out;
  //double dxsig_n_fac_out;
  //double dxsig_p_fac_out;
  //double nsigx_fid_out;
  //double nsigy_fid_out;
  //double W2low_out;
  //double W2high_out;
  //int passW2_out;
  //int passCoin_out;
  //int passFid_out;
  //setup new output tree branches
  //Parse->Branch("coin_mean", &coin_mean_out, "coin_mean/D");
  //Parse->Branch("coin_sigma", &coin_sigma_out, "coin_sigma/D");
  //Parse->Branch("dyO_p", &dyO_p_out, "dyO_p/D");
  //Parse->Branch("dyO_n", &dyO_n_out, "dyO_n/D");
  //Parse->Branch("dysig_p", &dysig_p_out, "dysig_p/D");
  //Parse->Branch("dysig_n", &dysig_n_out, "dysig_n/D");
  //Parse->Branch("dysig_n_fac", &dysig_n_fac_out, "dysig_n_fac/D");
  //Parse->Branch("dysig_p_fac", &dysig_p_fac_out, "dysig_p_fac/D");
  //Parse->Branch("dxO_p", &dxO_p_out, "dxO_p/D");
  //Parse->Branch("dxO_n", &dxO_n_out, "dxO_n/D");
  //Parse->Branch("dxsig_p", &dxsig_p_out, "dxsig_p/D");
  //Parse->Branch("dxsig_n", &dxsig_n_out, "dxsig_n/D");
  //Parse->Branch("dxsig_n_fac", &dxsig_n_fac_out, "dxsig_n_fac/D");
  //Parse->Branch("dxsig_p_fac", &dxsig_p_fac_out, "dxsig_p_fac/D");
  //Parse->Branch("nsigx_fid", &nsigx_fid_out , "nsigx_fid/D");
  //Parse->Branch("nsigy_fid", &nsigy_fid_out , "nsigy_fid/D");
  //Parse->Branch("W2low", &W2low_out, "W2low/D");
  //Parse->Branch("W2high", &W2high_out, "W2high/D");
  //Parse->Branch("passW2", &passW2_out,"passW2/I");
  //Parse->Branch("passCoin", &passCoin_out,"passCoin/I");
  //Parse->Branch("passFid", &passFid_out, "passFid/I");

  //loop over the run numbers
  for(int j = 0; j<num_runs; j++){

  //get run specific information
  data_object datData = myData[j];
  int run = datData.getRun();
  int field = datData.getSBSField();
  double Ebeam = datData.getBeamEnergy();
  double hcaldist = datData.getHCalDist();
  double hcaltheta = datData.getHCalAngle_Rad();
  double bbtheta = datData.getBBAngle_Rad();
  double sbsdist = datData.getSBSDist();
  TString input_file_name = datData.getInputFile();

  //add the file to the TChain
  C = new TChain("T");

  C->Add(input_file_name);

  // setting up ROOT tree branch addresses
  C->SetBranchStatus("*",0);

  //HCal general branches

  double sbs_hcal_x,sbs_hcal_y,sbs_hcal_e,sbs_hcal_nclus,sbs_hcal_index,sbs_hcal_nblk;

  C->SetBranchStatus("sbs.hcal.nblk",1);
  C->SetBranchStatus("sbs.hcal.x",1);
  C->SetBranchStatus("sbs.hcal.y",1);
  C->SetBranchStatus("sbs.hcal.e",1);
  C->SetBranchStatus("sbs.hcal.nclus",1);
  C->SetBranchStatus("sbs.hcal.index",1);

  C->SetBranchAddress("sbs.hcal.nblk",&sbs_hcal_nblk);
  C->SetBranchAddress("sbs.hcal.x", &sbs_hcal_x);
  C->SetBranchAddress("sbs.hcal.y", &sbs_hcal_y);
  C->SetBranchAddress("sbs.hcal.e", &sbs_hcal_e);
  C->SetBranchAddress("sbs.hcal.nclus", &sbs_hcal_nclus);
  C->SetBranchAddress("sbs.hcal.index", &sbs_hcal_index);

  //HCal cluster branches

  double sbs_hcal_clus_adctime[exp_constants::maxclus], sbs_hcal_clus_e[exp_constants::maxclus], sbs_hcal_clus_x[exp_constants::maxclus], sbs_hcal_clus_y[exp_constants::maxclus],sbs_hcal_clus_nblk[exp_constants::maxclus], sbs_hcal_clus_blk_row[exp_constants::maxclus], sbs_hcal_clus_blk_col[exp_constants::maxclus], sbs_hcal_clus_blk_id[exp_constants::maxclus], sbs_hcal_atimeblk;
  int Ndata_sbs_hcal_clus_id;

  C->SetBranchStatus("sbs.hcal.clus.adctime",1);
  C->SetBranchStatus("sbs.hcal.atimeblk",1);
  C->SetBranchStatus("sbs.hcal.clus.e",1);
  C->SetBranchStatus("sbs.hcal.clus.x",1);
  C->SetBranchStatus("sbs.hcal.clus.y",1);
  C->SetBranchStatus("sbs.hcal.clus.e",1);
  C->SetBranchStatus("Ndata.sbs.hcal.clus.id", 1);
  C->SetBranchStatus("sbs.hcal.clus.nblk", 1);
  C->SetBranchStatus("sbs.hcal.clus_blk.row", 1);
  C->SetBranchStatus("sbs.hcal.clus_blk.col", 1);
  C->SetBranchStatus("sbs.hcal.clus_blk.id", 1);

  C->SetBranchAddress("sbs.hcal.clus.adctime", &sbs_hcal_clus_adctime);
  C->SetBranchAddress("sbs.hcal.atimeblk", &sbs_hcal_atimeblk);
  C->SetBranchAddress("sbs.hcal.clus.e", &sbs_hcal_clus_e);
  C->SetBranchAddress("sbs.hcal.clus.x", &sbs_hcal_clus_x);
  C->SetBranchAddress("sbs.hcal.clus.y", &sbs_hcal_clus_y);
  C->SetBranchAddress("Ndata.sbs.hcal.clus.id", &Ndata_sbs_hcal_clus_id);
  C->SetBranchAddress("sbs.hcal.clus.nblk", &sbs_hcal_clus_nblk);
  C->SetBranchAddress("sbs.hcal.clus_blk.row", &sbs_hcal_clus_blk_row);
  C->SetBranchAddress("sbs.hcal.clus_blk.col", &sbs_hcal_clus_blk_col);
  C->SetBranchAddress("sbs.hcal.clus_blk.id", &sbs_hcal_clus_blk_id);

  //BBCal shower

  double bb_sh_atimeblk, bb_sh_e, bb_sh_nclus, bb_sh_nblk;

  C->SetBranchStatus( "bb.sh.atimeblk", 1 );
  C->SetBranchStatus( "bb.sh.e", 1 );
  C->SetBranchStatus( "bb.sh.nclus", 1 );
  C->SetBranchStatus( "bb.sh.nblk", 1 );

  C->SetBranchAddress("bb.sh.atimeblk", &bb_sh_atimeblk);
  C->SetBranchAddress("bb.sh.e", &bb_sh_e);
  C->SetBranchAddress("bb.sh.nclus", &bb_sh_nclus);
  C->SetBranchAddress("bb.sh.nblk", &bb_sh_nblk);

  //BBCal preshower

  double bb_ps_atimeblk, bb_ps_e, bb_ps_nclus, bb_ps_nblk;

  C->SetBranchStatus( "bb.ps.atimeblk", 1 );
  C->SetBranchStatus( "bb.ps.e", 1 );
  C->SetBranchStatus( "bb.ps.nclus", 1 );
  C->SetBranchStatus( "bb.ps.nblk", 1 );

  C->SetBranchAddress("bb.ps.atimeblk", &bb_ps_atimeblk);
  C->SetBranchAddress("bb.ps.e", &bb_ps_e);
  C->SetBranchAddress("bb.ps.nclus", &bb_ps_nclus);
  C->SetBranchAddress("bb.ps.nblk", &bb_ps_nblk);

  //BBGEM hits

  double bb_gem_track_nhits[maxtracks],bb_gem_track_ngoodhits[maxtracks],bb_gem_track_chi2ndf[maxtracks];

  C->SetBranchStatus("bb.gem.track.nhits", 1);
  C->SetBranchStatus("bb.gem.track.ngoodhits",1);
  C->SetBranchStatus("bb.gem.track.chi2ndf",1);

  C->SetBranchAddress("bb.gem.track.nhits",&bb_gem_track_nhits);
  C->SetBranchAddress("bb.gem.track.ngoodhits",&bb_gem_track_ngoodhits);
  C->SetBranchAddress("bb.gem.track.chi2ndf",&bb_gem_track_chi2ndf);

  // track branches

  double bb_tr_n, bb_tr_px[maxtracks], bb_tr_py[maxtracks], bb_tr_pz[maxtracks], bb_tr_p[maxtracks], bb_tr_x[maxtracks], bb_tr_y[maxtracks], bb_tr_vx[maxtracks], bb_tr_vy[maxtracks], bb_tr_vz[maxtracks], bb_tr_r_x[maxtracks], bb_tr_r_y[maxtracks],bb_tr_r_th[maxtracks], bb_tr_r_ph[maxtracks], bb_tr_th[maxtracks], bb_tr_ph[maxtracks];

  C->SetBranchStatus("bb.tr.n",1);
  C->SetBranchStatus("bb.tr.px",1);
  C->SetBranchStatus("bb.tr.py",1);
  C->SetBranchStatus("bb.tr.pz",1);
  C->SetBranchStatus("bb.tr.p",1);
  C->SetBranchStatus("bb.tr.x",1);
  C->SetBranchStatus("bb.tr.y",1);
  C->SetBranchStatus("bb.tr.vx",1);
  C->SetBranchStatus("bb.tr.vy",1);
  C->SetBranchStatus("bb.tr.vz",1);
  C->SetBranchStatus("bb.tr.th",1);
  C->SetBranchStatus("bb.tr.ph",1);
  C->SetBranchStatus("bb.tr.r_th",1);
  C->SetBranchStatus("bb.tr.r_ph",1);
  C->SetBranchStatus("bb.tr.r_x",1);
  C->SetBranchStatus("bb.tr.r_y",1);

  C->SetBranchAddress("bb.tr.n",&bb_tr_n);
  C->SetBranchAddress("bb.tr.px",&bb_tr_px);
  C->SetBranchAddress("bb.tr.py",&bb_tr_py);
  C->SetBranchAddress("bb.tr.pz",&bb_tr_pz);
  C->SetBranchAddress("bb.tr.p",&bb_tr_p);
  C->SetBranchAddress("bb.tr.x",&bb_tr_x);
  C->SetBranchAddress("bb.tr.y",&bb_tr_y);
  C->SetBranchAddress("bb.tr.vx",&bb_tr_vx);
  C->SetBranchAddress("bb.tr.vy",&bb_tr_vy);
  C->SetBranchAddress("bb.tr.vz",&bb_tr_vz);
  C->SetBranchAddress("bb.tr.th",&bb_tr_th);
  C->SetBranchAddress("bb.tr.ph",&bb_tr_ph);
  C->SetBranchAddress("bb.tr.r_th",&bb_tr_r_th);
  C->SetBranchAddress("bb.tr.r_ph",&bb_tr_r_ph);
  C->SetBranchAddress("bb.tr.r_x",&bb_tr_r_x);
  C->SetBranchAddress("bb.tr.r_y",&bb_tr_r_y);

  //grinch branches

  double bb_grinch_tdc_clus_trackindex, bb_grinch_tdc_clus_size, bb_grinch_tdc_hit_time[maxtracks], bb_grinch_tdc_hit_pmtnum[maxtracks], bb_grinch_tdc_hit_amp[maxtracks];
  int Ndata_bb_grinch_tdc_hit_time;

  C->SetBranchStatus("bb.grinch_tdc.clus.trackindex", 1);
  C->SetBranchStatus("bb.grinch_tdc.clus.size", 1);
  C->SetBranchStatus("bb.grinch_tdc.hit.time", 1);
  C->SetBranchStatus("bb.grinch_tdc.hit.pmtnum", 1);
  C->SetBranchStatus("bb.grinch_tdc.hit.amp", 1);
  C->SetBranchStatus("Ndata.bb.grinch_tdc.hit.time", 1);

  C->SetBranchAddress("bb.grinch_tdc.clus.trackindex", &bb_grinch_tdc_clus_trackindex);
  C->SetBranchAddress("bb.grinch_tdc.clus.size", &bb_grinch_tdc_clus_size);
  C->SetBranchAddress("bb.grinch_tdc.hit.time", &bb_grinch_tdc_hit_time);
  C->SetBranchAddress("bb.grinch_tdc.hit.pmtnum", &bb_grinch_tdc_hit_pmtnum);
  C->SetBranchAddress("bb.grinch_tdc.hit.amp", &bb_grinch_tdc_hit_amp);
  C->SetBranchAddress("Ndata.bb.grinch_tdc.hit.time", &Ndata_bb_grinch_tdc_hit_time);

  //ekine branches

  double e_kine_Q2, e_kine_W2, e_kine_epsilon, e_kine_nu, e_kine_q_x, e_kine_q_y, e_kine_q_z;

  C->SetBranchStatus("e.kine.Q2",1);
  C->SetBranchStatus("e.kine.W2",1);
  C->SetBranchStatus("e.kine.epsilon",1);
  C->SetBranchStatus("e.kine.nu",1);
  C->SetBranchStatus("e.kine.q_x",1);
  C->SetBranchStatus("e.kine.q_y",1);
  C->SetBranchStatus("e.kine.q_z",1);

  C->SetBranchAddress("e.kine.Q2", &e_kine_Q2);
  C->SetBranchAddress("e.kine.W2", &e_kine_W2);
  C->SetBranchAddress("e.kine.epsilon", &e_kine_epsilon);
  C->SetBranchAddress("e.kine.nu", &e_kine_nu);
  C->SetBranchAddress("e.kine.q_x", &e_kine_q_x);
  C->SetBranchAddress("e.kine.q_y", &e_kine_q_y);
  C->SetBranchAddress("e.kine.q_z", &e_kine_q_z);

  //CODA event variables

  double g_trigbits, g_evtime;

  C->SetBranchStatus("g.trigbits",1);
  C->SetBranchStatus("g.evtime",1);

  C->SetBranchAddress("g.trigbits", &g_trigbits);
  C->SetBranchAddress("g.evtime", &g_evtime);

  //global cut branches
  //already handled above

  //setup global cut formula
  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", globalcut, C );

  //setup hcal coordinate system with hcal angle with respect to exist beamline
  TVector3 hcal_zaxis = physics::getHCal_zaxis(hcaltheta);
  TVector3 hcal_xaxis = physics::getHCal_xaxis();
  TVector3 hcal_yaxis = physics::getHCal_yaxis(hcal_xaxis,hcal_zaxis);
  //define HCal origin
  TVector3 hcal_origin = physics::getHCal_origin(hcaldist,hcal_offset,hcal_xaxis,hcal_zaxis);
  //mean energy loss of the beam before scattering
  double Eloss_outgoing = physics::getEloss_outgoing(bbtheta,target);

  //Accounting event variables
  long nevent = 0, nentries = C->GetEntries();

  //ttree formula variables
  int treenum = 0, currenttreenum = 0;

  	//event loop
  	while(C->GetEntry(nevent++)){

	//progress tracker
	cout << "Processing run " <<  j << "/" << num_runs << " run number " << run << " event " << nevent << "/" << nentries << "\r";
	cout.flush();

	//single loop global cut
	currenttreenum = C->GetTreeNumber();
    	if( nevent == 1 || currenttreenum != treenum ){
    	  treenum = currenttreenum;
    	  GlobalCut->UpdateFormulaLeaves();

          auto* Run_Data = C->GetFile()->Get<THaRunBase>("Run_Data");
          TDatime run_time = Run_Data->GetDate();
          run_time.Set(run_time.GetYear(),run_time.GetMonth(),run_time.GetDay(),run_time.GetHour(),run_time.GetMinute(),0);
          run_time_unix = run_time.Convert();
    	}
        //Is true if failed global cut
    	bool failglobal = cuts::failedGlobal(GlobalCut);

        double time_interval = 4; //in ns -> shouldn't this be 2 for GEn?
        int time_rel = g_evtime*time_interval*1e-9/60; // in min, rounded
	TDatime time_abs(run_time_unix + time_rel * 60);

	auto it = DBInfo.He3Pol.find(time_abs);
	if(it == DBInfo.He3Pol.end())
	  He3Pol_out = -1;
	else
	  He3Pol_out = it->second;

	datetime_out = time_abs;

	///////////
	//Electron-arm physics calculations

	//BB E/p calculated. Should be the same as the tree variable. But it does not cause seg faults
	double BB_E_over_p = (bb_sh_e+bb_ps_e)/bb_tr_p[0];

	//correct beam energy from vertex information
	double Eloss = physics::getEloss(bb_tr_vz[0],target);
	double Ecorr = physics::getEcorr(Ebeam,Eloss);

	//make the vertex, assuming only beam direction
	TVector3 vertex = physics::getVertex(bb_tr_vz[0]);

	//reconstructed momentum from track momentum information, corrected for mean energy loss. Still need to include losses from Al shielding or target windows later
	double pcorr = physics::getp_recon_corr(bb_tr_p[0],Eloss_outgoing);

        //four momentum vector for electron beam with correted Energy value
        TLorentzVector pbeam = physics::getpBeam(Ecorr);

	//four momentum for scattered electron based on reconstruction
	TLorentzVector p_eprime = physics::getp_eprime(bb_tr_px[0],bb_tr_py[0],bb_tr_pz[0],bb_tr_p[0],pcorr);

	//four vector for target
	TLorentzVector p_targ = physics::getp_targ(target);

	//four vector, virtual photon momentum or momentum transferred to the scattered nucleon
	TLorentzVector q = physics::getq(pbeam,p_eprime);
	TVector3 q_vec = q.Vect();

	//Theta for scattered electron using reconstructed track momentum
	double etheta = physics::get_etheta(p_eprime);

	//Phi for scattered electron using reconstructed track momentum
	double ephi = physics::get_ephi(p_eprime);

	//central momentum reconstructed from track angles and beam energy
	double pcentral = physics::get_pcentral(pbeam,etheta,target);

	//assume coplanarity, get the expected phi for the nucleon
	double phi_N_exp = physics::get_phinucleon(ephi,physics_constants::PI);

	//Calculate Mott cross section for this event
	double Mott_CS = physics::getMott_CS(physics_constants::alpha,etheta,pcorr,Ecorr);

	/* Can reconstruct e' momentum for downstream calculations differently:
	* v1 - Use four-momentum member functions
 	* v2 - Use all available ekine (tree) vars and calculate vectors (should be the same as v1)
	* v3 - Use reconstructed angles as independent qty (usually preferable given GEM precision at most kinematics)
 	* v4 - Use reconstructed momentum as independent qty */


	//four momentum transferred squared
	double Q2;

  	//four vector, scattered nucleon momentum
  	TLorentzVector p_N;
	TVector3 p_Nhat;

	//scattered nucleon expected momentum
	double p_N_exp;

	//energy transfer
	double nu;

	//scattered nucleon expected angle theta
	double theta_N_exp;

	//Invariant Mass Squared
	double W2;

	//scaling variable tau
	double tau;

	//polarization of the virtual photon
	double epsilon;

	//conditional to determine remaining e-arm related calculations.
		if(e_method == 1){
		//v1
		Q2 = physics::getQ2(q);
		p_N = physics::get_pN(q,p_targ);
		p_Nhat = p_N.Vect().Unit();
		nu = physics::getnu(q);
		W2 = physics::getW2(p_N);
		tau = physics::get_tau(Q2,target);
		epsilon = physics::get_epsilon(tau,etheta);
		}else if(e_method == 2){
		//v2
		Q2 = physics::getQ2(e_kine_Q2);
		W2 = physics::getW2(e_kine_W2);
		nu = physics::getnu(e_kine_nu);
		p_N_exp = physics::get_pNexp(nu,target);
		theta_N_exp = physics::get_thetaNexp(pbeam,pcentral,etheta,p_N_exp);
		p_Nhat = physics::get_pNhat(theta_N_exp,phi_N_exp);
		p_N = physics::get_pN(p_N_exp,p_Nhat,nu,p_targ);
		tau = physics::get_tau(Q2,target);
		epsilon = physics::get_epsilon(tau,etheta);
		}else if(e_method == 3){
		//v3
		Q2 = physics::getQ2(pbeam,p_eprime,etheta);
		nu = physics::getnu(pbeam,pcentral);
		p_N_exp = physics::get_pNexp(nu,target);
                theta_N_exp = physics::get_thetaNexp(pbeam,pcentral,etheta,p_N_exp);
		p_Nhat = physics::get_pNhat(theta_N_exp,phi_N_exp);
                p_N = physics::get_pN(p_N_exp,p_Nhat,nu,p_targ);
		W2 = physics::getW2(pbeam,p_eprime,Q2,target);
		tau = physics::get_tau(Q2,target);
		epsilon = physics::get_epsilon(tau,etheta);
		}else if(e_method == 4){
		//v4
		Q2 = physics::getQ2(pbeam,p_eprime,etheta);
		nu = physics::getnu(pbeam,p_eprime);
		p_N_exp = physics::get_pNexp(nu,target);
                theta_N_exp = physics::get_thetaNexp(pbeam,pcentral,etheta,p_N_exp);
                p_Nhat = physics::get_pNhat(theta_N_exp,phi_N_exp);
                p_N = physics::get_pN(p_N_exp,p_Nhat,nu,p_targ);
                W2 = physics::getW2(pbeam,p_eprime,Q2,target);
		tau = physics::get_tau(Q2,target);
		epsilon = physics::get_epsilon(tau,etheta);
		}else{
		//Error handling, default version 3
		cout << "Warning: Method for calculating e-arm physics was not included. Defaulting to method 3." << endl;
		Q2 = physics::getQ2(pbeam,p_eprime,etheta);
		nu = physics::getnu(pbeam,pcentral);
                p_N_exp = physics::get_pNexp(nu,target);
                theta_N_exp = physics::get_thetaNexp(pbeam,pcentral,etheta,p_N_exp);
                p_Nhat = physics::get_pNhat(theta_N_exp,phi_N_exp);
                p_N = physics::get_pN(p_N_exp,p_Nhat,nu,p_targ);
                W2 = physics::getW2(pbeam,p_eprime,Q2,target);
		tau = physics::get_tau(Q2,target);
		epsilon = physics::get_epsilon(tau,etheta);
		}

	// ray from Hall origin onto the face of hcal where the nucleon hit
	TVector3 hcal_intersect = physics::get_hcalintersect(vertex,hcal_origin,hcal_zaxis,p_Nhat );

	//gets expected location of scattered nucleon assuming straight line projections from BB track, x-direction
	double xhcal_expect = physics::get_xhcalexpect(hcal_intersect,hcal_origin,hcal_xaxis);

	//gets expected location of scattered nucleon assuming straight line projections from BB track, y-direction
        double yhcal_expect = physics::get_yhcalexpect(hcal_intersect,hcal_origin,hcal_yaxis);

	//Calculate expected proton deflection with a somewhat crude module
	double BdL = physics::getBdL(sbs_field,kin,pass);
	double proton_deflection = physics::get_protonDeflection(BdL,p_N.Vect().Mag(),hcaldist,sbsdist);

	//////////////////////
	//INTIME CLUSTER ANALYSIS
	//Requires that it has the greatest hcal cluster energy and that hcal cluster analog time is in coincidence with bbcal analog time

	//intime cluster selection analysis, intime algorithm
	//int intime_idx = physics::cluster_intime_select(Ndata_sbs_hcal_clus_id,sbs_hcal_clus_adctime,bb_sh_atimeblk,sbs_hcal_clus_e,coin_mean,coin_sig_fac,coin_profile_sig,hcalemin);

	//Assume that the itime analysis is sufficient to find the best cluster in HCal
	//int clus_idx_best = intime_idx;

	//calculate important information from best cluster
  // HCAL FIRST CLUSTER IS ALREADY "BEST" CLUSTER, NO NEED TO DEFINE IT
	//double xhcal_bestclus = sbs_hcal_clus_x[clus_idx_best];
 	//double yhcal_bestclus = sbs_hcal_clus_y[clus_idx_best];
	//double xhcal_pclus = sbs_hcal_clus_x[0];
  //double yhcal_pclus = sbs_hcal_clus_y[0];
	//double dx_bestclus = physics::get_dx(xhcal_bestclus,xhcal_expect);
	//double dx_pclus = physics::get_dx(xhcal_pclus,xhcal_expect);
 	//double dy_bestclus = physics::get_dy(yhcal_bestclus,yhcal_expect);

  double dx = physics::get_dx(sbs_hcal_x,xhcal_expect);
  double dy = physics::get_dx(sbs_hcal_y,yhcal_expect);

	//double hcal_atime_bestclus = sbs_hcal_clus_adctime[clus_idx_best];
	double coin = sbs_hcal_atimeblk - bb_sh_atimeblk;
	//double coin_pclus = sbs_hcal_clus_adctime[0] - bb_sh_atimeblk;
	//double hcal_e_bestclus = sbs_hcal_clus_e[clus_idx_best];
	//int hcal_nblk_bestclus = (int) sbs_hcal_clus_nblk[clus_idx_best];

	//calculate the number of sigma away from fiducial boundaries. Store info for later
        //double nsigx_fid = cuts::calculate_nsigma_fid_x(xhcal_expect,dxsig_p,dxsig_n,dx_pn,hcalaa);
        //double nsigy_fid = cuts::calculate_nsigma_fid_y(yhcal_expect,dysig_p,hcalaa);


	//setup booleans for cuts later. Save boolean values to tree
	//global is above

	//HCal active area
	bool hcalaa_ON = cuts::hcalaa_ON(sbs_hcal_x,sbs_hcal_y,hcalaa);
	bool hcalaa_ON_exp = cuts::hcalaa_ON(xhcal_expect,yhcal_expect,hcalaa);

	//W2 elastic boolean
	bool goodW2 = cuts::goodW2(W2,W2_low,W2_high);
	bool antiW2_low = W2 > 0.0 && W2 < W2_low;
	bool antiW2_high = W2 > W2_high;
	bool antiW2_morehigh = W2 > (W2_high+1.1);

	//good dy boolean
	//bool good_dy = cuts::good_dy(dy_bestclus,dyO_p,dysig_cut_fac,dysig_cut);

	//good coincidence time cut
	//bool passCoin = cuts::passCoin(coin_bestclus,coin_mean,coin_sig_fac,coin_sig);
	//Used for coin anti cut
	//bool passCoin_pclus = cuts::passCoin(coin_pclus,coin_mean,coin_sig_fac,coin_sig);

	//good fiducial cut
	//bool passFid = cuts::hcalfid_IN(xhcal_expect,yhcal_expect,dx_pn,hcalfid);

	//pass HCal E
	//bool passHCalE = cuts::passHCalE(hcal_e_bestclus,hcalemin);

	//pass HCal num clus
	//bool passHCal_Nclus = cuts::passHCal_NClus(sbs_hcal_nclus,hcalnclusmin);

	//pass NSig Fid check
	//bool passNSigFid = cuts::passNsigFid(nsigx_fid,nsigy_fid);

	//Fill analysis tree variables before making cuts
	//dx_out = dx_bestclus;
  dx_out = dx;
 	//dy_out = dy_bestclus;
  dy_out = dy;
  xexp_out = xhcal_expect;
  yexp_out = yhcal_expect;
  sbs_hcal_x_out = sbs_hcal_x;
 	sbs_hcal_y_out = sbs_hcal_y;
  W2_out = W2;
  nu_out = nu;
	tau_out = tau;
	epsilon_out = epsilon;
	pcorr_out = pcorr;
	mott_out = Mott_CS;

  sbs_hcal_e_out = sbs_hcal_e;

  BBtot_e_out = bb_sh_e+bb_ps_e;
  bb_sh_e_out = bb_sh_e;
  bb_ps_e_out = bb_ps_e;

  sbs_hcal_atimeblk_out = sbs_hcal_atimeblk;

  bb_sh_atimeblk_out = bb_sh_atimeblk;
  bb_ps_atimeblk_out = bb_ps_atimeblk;

  bb_gem_track_nhits_out = bb_gem_track_nhits[0];
  bb_gem_track_ngoodhits_out = bb_gem_track_ngoodhits[0];
	bb_gem_track_chi2ndf_out = bb_gem_track_chi2ndf[0];

  bb_tr_x_out = bb_tr_x[0];
  bb_tr_y_out = bb_tr_y[0];
  bb_tr_p_out = bb_tr_p[0];
  bb_tr_vz_out = bb_tr_vz[0];

  BB_E_over_p_out = BB_E_over_p;

  bb_tr_th_out = bb_tr_th[0];
	bb_tr_ph_out = bb_tr_ph[0];
	bb_tr_r_x_out = bb_tr_r_x[0];
	bb_tr_r_y_out = bb_tr_r_y[0];
	bb_tr_r_th_out = bb_tr_r_th[0];
	bb_tr_r_ph_out = bb_tr_r_ph[0];

	Ndata_sbs_hcal_clus_id_out = Ndata_sbs_hcal_clus_id ;
  sbs_hcal_nclus_out = sbs_hcal_nclus;

  Ndata_bb_grinch_tdc_hit_time_out = Ndata_bb_grinch_tdc_hit_time;
  bb_grinch_tdc_clus_size_out = bb_grinch_tdc_clus_size;
  bb_grinch_tdc_clus_trackindex_out = bb_grinch_tdc_clus_trackindex;
  for(int iclus = 0; iclus < Ndata_bb_grinch_tdc_hit_time; iclus++)
  {
    bb_grinch_tdc_hit_time_out[iclus] = bb_grinch_tdc_hit_time[iclus];
    bb_grinch_tdc_hit_pmtnum_out[iclus] = bb_grinch_tdc_hit_pmtnum[iclus];
    bb_grinch_tdc_hit_amp_out[iclus] = bb_grinch_tdc_hit_amp[iclus];
  }

  bb_sh_nclus_out = (int) bb_sh_nclus;
	bb_sh_nblk_out = (int) bb_sh_nblk;
	bb_ps_nclus_out = (int) bb_ps_nclus;
  bb_ps_nblk_out = (int) bb_ps_nblk;
	bb_tr_n_out = bb_tr_n;

  passGlobal_out = (int) !failglobal;
	HCalON_out = (int) hcalaa_ON;

  g_evtime_out = g_evtime;
  g_trigbits_out = g_trigbits;

  //date time???

  run_out = run ;
  mag_out = field;

  hcal_sh_atime_diff_out = coin;
	proton_deflection_out = proton_deflection;
	p_central_out = pcentral;

  sbs_hcal_nblk_out = sbs_hcal_nblk;
	//sbs_hcal_clus_blk_row_out = sbs_hcal_clus_blk_row[clus_idx_best];
  sbs_hcal_clus_blk_row_out = sbs_hcal_clus_blk_row[0];
	//sbs_hcal_clus_blk_col_out = sbs_hcal_clus_blk_col[clus_idx_best];
  sbs_hcal_clus_blk_col_out = sbs_hcal_clus_blk_col[0];
  sbs_hcal_clus_blk_id_out = (int) sbs_hcal_clus_blk_id[0];

	//For physics extraction
	Ebeam_corr_out = Ecorr;
	E_eprime_out =  p_eprime.E();
	etheta_out = etheta;
	p_N_out = p_N.Vect().Mag();
	Q2_out = Q2;

	//Fill the analysis tree
	Parse->Fill();

	}//end event loop

	// reset chain for the next run config
	C->Reset();

  }//end loop over the run numbers

  //Write everything to output file
  fout->Write();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << watch->CpuTime() << " s = " << watch->CpuTime()/60.0 << " min. Real time = " << watch->RealTime() << " s = " << watch->RealTime()/60.0 << " min." << endl;

}//end Main
