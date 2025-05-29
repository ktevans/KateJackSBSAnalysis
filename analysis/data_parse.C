//Author Ezekiel Wertz
//02/01/2024
//Purpose: Parsing Script for LD2 data to produce output histograms for later analysis

//The exact ordering of this matters. ROOT for some reason cannot handle calls for files that have already been included. 
#include "TF1.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCut.h"
#include "TEventList.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "../../src/utility.C"
#include "../../src/exp_constants.C"
#include "../../src/kinematic_obj.C"
#include "../../src/fits.C"
#include "../../src/data_object.C"
#include "../../src/cuts.C"
#include "../../src/physics.C"
#include "../../src/parse_config.C"
#include "../../src/plots.C"
#include "../../src/calc_FFs_RCS_obj.C"
#include "../../src/fit_histogram.C"

//Main
void data_elastic_parse(const char *setup_file_name){

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
  double hcalfit_low = exp_constants::hcalposXi_mc; //lower fit/bin limit for hcal dx plots. 
  double hcalfit_high = exp_constants::hcalposXf_mc; //higher fit/bin limit for hcal dx plots.
  double hcal_fitrange = exp_constants::hcal_vrange; //Full range of hcal dx plots

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

  //Histograms///////

  //global cuts
  TH1D *h_ntracks = new TH1D("ntracks","Number of Tracks;", 150, 0, 5);
  TH1D *h_ntracks_globcut = new TH1D("ntracks_globcut","Number of Tracks,global cut;", 150, 0, 5);
  TH1D *h_ntracks_cut = new TH1D("ntracks_cut","Number of Tracks, cuts;", 150, 0, 5);
  TH1D *h_track_chi2ndf_nocut = new TH1D("track_chi2ndf_nocut","Track Chi^2/NDF, no cuts;",200,0,50);
  TH1D *h_track_chi2ndf_globcut = new TH1D("track_chi2ndf_globcut","Track Chi^2/NDF, global cuts;",200,0,20);
  TH1D *h_track_chi2ndf_cut = new TH1D("track_chi2ndf_cut","Track Chi^2/NDF, all cuts;",200,0,20);
  TH1D *h_PS_E = new TH1D("h_ps_e"," PS Cluster Energy (GeV);",250,0.0,2.5);
  TH1D *h_PS_E_globcut = new TH1D("h_ps_e_globcut"," PS Cluster Energy (GeV),global cut;",250,0.0,2.5);
  TH1D *h_PS_E_cut = new TH1D("h_ps_e_cut"," PS Cluster Energy (GeV), Cuts;",250,0.0,2.5);
  TH1D *h_vert_z = new TH1D( "vert_z", "Vertex Position z-direction (m); m", 200, -0.15, 0.15 );
  TH1D *h_vert_z_globcut = new TH1D( "vert_z_globcut", "Vertex Position z-direction (m), global cut; m", 200, -0.15, 0.15 );
  TH1D *h_vert_z_cut = new TH1D( "vert_z_cut", "Vertex Position z-direction (m), Cuts; m", 200, -0.15, 0.15 );
  TH1D *h_HCal_E = new TH1D( "HCal_E", "HCal Cluster Energy (GeV); GeV", 250, 0, 0.4 );
  TH1D *h_HCal_E_globcut = new TH1D( "HCal_E_globcut", "HCal Cluster Energy (GeV), global cut; GeV", 250, 0, 0.4 );
  TH1D *h_HCal_E_cut = new TH1D( "HCal_E_cut", "HCal Cluster Energy (GeV), Cuts; GeV", 500, 0, 1.0 );
  TH1D *h_HCal_E_cut_noHCalE = new TH1D( "HCal_E_cut_noHCalE", "HCal Cluster Energy (GeV), all cuts but HCal E; GeV", 500, 0, 1.0 );
  TH1D *h_HCal_E_best = new TH1D( "HCal_E_best", "HCal Cluster Energy (GeV), best cluster; GeV", 250, 0, 0.4 );
  TH1D *h_HCal_E_best_globcut = new TH1D( "HCal_E_best_globcut", "HCal Cluster Energy (GeV), best cluster, global cut; GeV", 250, 0, 0.4 );
  TH1D *h_HCal_E_best_cut = new TH1D( "HCal_E_best_cut", "HCal Cluster Energy (GeV), best cluster, Cuts; GeV", 500, 0, 1.0 );
  TH1D *h_HCal_E_best_cut_noHCalE = new TH1D( "HCal_E_best_cut_noHCalE", "HCal Cluster Energy (GeV), best cluster,all  Cuts but HCal E; GeV", 500, 0, 1.0 );
  TH1D *h_HCal_nclus = new TH1D("HCal_nclus","HCal number of clusters meeting threshold;", 250,0,10);
  TH1D *h_HCal_nclus_globcut = new TH1D("HCal_nclus_globcut","HCal number of clusters meeting threshold, global cut;", 250,0,10);
  TH1D *h_HCal_nclus_cut = new TH1D("HCal_nclus_cut","HCal number of clusters meeting threshold, Cuts;", 250,0,10);
  TH1D *h_SH_nclus = new TH1D("SH_nclus","SH number of clusters meeting threshold;", 250,0,10);
  TH1D *h_SH_nclus_globcut = new TH1D("SH_nclus_globcut","SH number of clusters meeting threshold, global cut;", 250,0,10);
  TH1D *h_SH_nclus_cut = new TH1D("SH_nclus_cut","SH number of clusters meeting threshold, Cuts;", 250,0,10);
  TH1D *h_TPS_SH_globcut = new TH1D("h_tps_sh_globcut","Total PS and SH cluster energy (GeV),global cut;",250,0.5,4.0);
  TH1D *h_TPS_SH_cut = new TH1D("h_tps_sh_cut","Total PS and SH cluster energy (GeV), All cuts;",250,0.5,4.0);
  TH1D *h_nhits = new TH1D("nhits","Number of hits on track;",150, 0, 6);
  TH1D *h_nhits_globcut = new TH1D("nhits_globcut","Number of hits on track, global cut;",150, 0, 6);
  TH1D *h_nhits_cut = new TH1D("nhits_cut","Number of hits on track, all cuts;",150, 0, 6);
  TH1D *h_bbtrp_nocut = new TH1D("bbtrp_nocut","BigBite Track Momentum (GeV), no cut;",300, 0.0, 4.0);
  TH1D *h_bbtrp_globcut = new TH1D("bbtrp_globcut","BigBite Track Momentum (GeV), global cut;",300, 0.0, 4.0);
  TH1D *h_bbtrp_cut = new TH1D("bbtrp_cut","BigBite Track Momentum (GeV), cuts;",300, 0.0, 4.0);
  TH1D *h_bbEoverp_nocut = new TH1D("bbEoverp_nocut","BigBite E over p, no cut;",300, 0.2, 1.8);
  TH1D *h_bbEoverp_globcut = new TH1D("bbEoverp_globcut","BigBite E over p, global cut;",300, 0.2, 1.8);
  TH1D *h_bbEoverp_cut = new TH1D("bbEoverp_cut","BigBite E over p, cuts;",300, 0.2, 1.8);
  TH1D *h_optics_xdir_globcut = new TH1D("h_optics_xdir_globcut", "BigBite optics validity, track x-dir;",200,-0.6,0.6);
  TH1D *h_optics_ydir_globcut = new TH1D("h_optics_ydir_globcut", "BigBite optics validity, track y-dir;",200,-0.2,0.2);
  TH2D *h_W2_optics_xdir_globcut = new TH2D("h_W2_optics_xdir_globcut", "W2 vs BigBite optics validity, track x-dir;", 200, -0.6, 0.6,binfac*W2fitmax, 0.0, W2fitmax );
  TH2D *h_W2_optics_ydir_globcut = new TH2D("h_W2_optics_ydir_globcut", "W2 vs BigBite optics validity, track y-dir;", 200, -0.2, 0.2,binfac*W2fitmax, 0.0, W2fitmax );
  TH2D *h_xexp_optics_xdir_globcut = new TH2D("h_xexp_optics_xdir_globcut", "HCal X Expect vs BigBite optics validity, track x-dir;", 200, -0.6, 0.6,600, -3.0, 3.0 );
  TH2D *h_yexp_optics_ydir_globcut = new TH2D("h_yexp_optics_ydir_globcut", "HCal Y Expect vs BigBite optics validity, track y-dir;", 200, -0.2, 0.2,400, -2.0, 2.0 );
  TH1D *h_optics_xdir_cut = new TH1D("h_optics_xdir_cut", "BigBite optics validity, track x-dir;",200,-0.6,0.6);
  TH1D *h_optics_ydir_cut = new TH1D("h_optics_ydir_cut", "BigBite optics validity, track y-dir;",200,-0.2,0.2);
  TH2D *h_W2_optics_xdir_cut = new TH2D("h_W2_optics_xdir_cut", "W2 vs BigBite optics validity, track x-dir;", 200, -0.6, 0.6,binfac*W2fitmax, 0.0, W2fitmax );
  TH2D *h_W2_optics_ydir_cut = new TH2D("h_W2_optics_ydir_cut", "W2 vs BigBite optics validity, track y-dir;", 200, -0.2, 0.2,binfac*W2fitmax, 0.0, W2fitmax );
  TH2D *h_xexp_optics_xdir_cut = new TH2D("h_xexp_optics_xdir_cut", "HCal X Expect vs BigBite optics validity, track x-dir;", 200, -0.6, 0.6,600, -3.0, 3.0 );
  TH2D *h_yexp_optics_ydir_cut = new TH2D("h_yexp_optics_ydir_cut", "HCal Y Expect vs BigBite optics validity, track y-dir;", 200, -0.2, 0.2,400, -2.0, 2.0 );
  TH1D *h_optics_xdir_nocut = new TH1D("h_optics_xdir_nocut", "BigBite optics validity, track x-dir;",200,-0.6,0.6);
  TH1D *h_optics_ydir_nocut = new TH1D("h_optics_ydir_nocut", "BigBite optics validity, track y-dir;",200,-0.2,0.2);
  TH2D *h_W2_optics_xdir_nocut = new TH2D("h_W2_optics_xdir_nocut", "W2 vs BigBite optics validity, track x-dir;", 200, -0.6, 0.6,binfac*W2fitmax, 0.0, W2fitmax );
  TH2D *h_W2_optics_ydir_nocut = new TH2D("h_W2_optics_ydir_nocut", "W2 vs BigBite optics validity, track y-dir;", 200, -0.2, 0.2,binfac*W2fitmax, 0.0, W2fitmax );
  TH2D *h_xexp_optics_xdir_nocut = new TH2D("h_xexp_optics_xdir_nocut", "HCal X Expect vs BigBite optics validity, track x-dir;", 200, -0.6, 0.6,600, -3.0, 3.0 );
  TH2D *h_yexp_optics_ydir_nocut = new TH2D("h_yexp_optics_ydir_nocut", "HCal Y Expect vs BigBite optics validity, track y-dir;", 200, -0.2, 0.2,400, -2.0, 2.0 );

  //basic H-arm
  TH2D *hxy_globcut = new TH2D("hxy_globcut","HCal X  vs Y, global cut;HCal Y  (m); HCal X  (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_glob_W2_cut = new TH2D("hxy_glob_W2_cut","HCal X  vs Y, global & W2 cuts;HCal Y  (m); HCal X  (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_acceptancecut = new TH2D("hxy_accpetancecut","HCal X  vs Y,global & W2 & acceptance cuts, ;HCal Y  (m); HCal X  (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_cut = new TH2D("hxy_cut","HCal X vs Y, all cuts, best cluster;y_{HCAL} (m); x_{HCAL} (m)",300, -2.0, 2.0, 500, -3.0, 3.0);

  //E-arm
  TH1D *h_W2_nocut = new TH1D( "W2_nocut", "W2 (GeV^{2}) no cut; GeV^{2}", binfac*5.0, 0.0, 5.0 );
  TH1D *h_W2_globcut = new TH1D( "W2_globcut", "W2 (GeV^{2}) global cut; GeV^{2}", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *h_W2_glob_W2_cut = new TH1D( "W2_glob_W2_cut", "W2 (GeV^{2}) global & W2 cuts; GeV^{2}", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *h_W2_cut = new TH1D( "W2_cut", "W2 (GeV^{2}) all cuts; GeV^{2}", binfac*W2fitmax, 0.0, W2fitmax );  
  TH1D *h_W2_notW2_cut = new TH1D( "W2_notW2_cut", "W2 (GeV^{2}) all cuts, but W2; GeV^{2}", binfac*5.0, 0.0, 5.0 );
  TH1D *h_Q2_globcut = new TH1D( "Q2_globcut", "Q2 (GeV^{2}) global cut; GeV^{2}", 300, 0.0, 6.0 );
  TH1D *hx_expect_nocut = new TH1D("hx_expect_nocut","HCal X Expect, no cut; HCal X Expect (m)", 600, -3.0, 3.0 );
  TH1D *hx_expect_globcut = new TH1D("hx_expect_globcut","HCal X Expect, global cut; HCal X Expect (m)", 600, -3.0, 3.0 );
  TH1D *hx_expect_cut = new TH1D("hx_expect_cut","HCal X Expect, all cuts; HCal X Expect (m)", 600, -3.0, 3.0 );
  TH1D *hy_expect_nocut = new TH1D("hy_expect_nocut","HCal Y Expect, no cut; HCal Y Expect (m)", 400, -2.0, 2.0 );
  TH1D *hy_expect_globcut = new TH1D("hy_expect_globcut","HCal Y Expect, global cut; HCal Y Expect (m)", 400, -2.0, 2.0 );
  TH1D *hy_expect_cut = new TH1D("hy_expect_cut","HCal Y Expect, all cuts; HCal Y Expect (m)", 400, -2.0, 2.0 );
  TH2D *hxy_expect_nocut = new TH2D("hxy_expect_nocut","HCal X Expect vs Y Expect, no cut;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_expect_nocut_n = new TH2D("hxy_expect_nocut_n","HCal X Expect vs Y Expect, no cut neu;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_expect_nocut_p = new TH2D("hxy_expect_nocut_p","HCal X Expect vs Y Expect, no cut pro;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_expect_globcut = new TH2D("hxy_expect_globcut","HCal X Expect vs Y Expect, global cut;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_expect_glob_W2_cut = new TH2D("hxy_expect_glob_W2_cut","HCal X Expect vs Y Expect, global & W2 cuts;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  
  TH2D *hxy_expect_glob_W2_cut_n = new TH2D("hxy_expect_glob_W2_cut_n","HCal X Expect vs Y Expect, global & W2 cuts neu;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );

  TH2D *hxy_expect_glob_W2_cut_p = new TH2D("hxy_expect_glob_W2_cut_p","HCal X Expect vs Y Expect, global & W2 cuts pro;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_expect_n = new TH2D("hxy_expect_n","HCal X Expect vs Y Expect, elastic cuts neutron;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_expect_p = new TH2D("hxy_expect_p","HCal X Expect vs Y Expect, elastic cuts proton;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_expect_fidcutp = new TH2D("hxy_expect_fidcutp","HCal X Expect vs Y Expect, proton passed fiducial;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
 TH2D *hxy_expect_fidcutn = new TH2D("hxy_expect_fidcutn","HCal X Expect vs Y Expect, neutron passed fiducial;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_expect_failedfid = new TH2D("hxy_expect_failedfid","HCal X Expect vs Y Expect, failed fiducial cut;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_expect_cut = new TH2D("hxy_expect_cut","HCal X Expect vs Y Expect, all cuts;HCal Y Expect (m); HCal X Expect (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );

  //Both arms
  TH2D *hdxdy_globcut = new TH2D("dxdy_globcut","HCal dxdy, global cut ;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",400,dy_low,dy_high,600,dx_low,dx_high);
  TH2D *hdxdy_glob_W2_cut = new TH2D("dxdy_glob_W2_cut","HCal dxdy, global & W2 cuts ;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",400,dy_low,dy_high,600,dx_low,dx_high);
  TH2D *hdxdy_cut = new TH2D("hdxdy_cut","HCal dxdy, all cuts;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",400,dy_low,dy_high,600,dx_low,dx_high);
  TH2D *hdxdy_nofid = new TH2D("dxdy_nofid","HCal dxdy, all cuts but fiducial ;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",400,dy_low,dy_high,600,dx_low,dx_high);
  TH2D *hdxdy_nody = new TH2D("hdxdy_nody","HCal dxdy, all cuts but dy ;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)",400,dy_low,dy_high,600,dx_low,dx_high);
  TH2D *hdxvE = new TH2D("dxvE","dx vs HCal E, all cuts;E_{HCAL} (GeV);x_{HCAL}-x_{expect} (m)", 400, 0.0, 4.0, 600, dx_low, dx_high );
  TH2D *hdxvW2 = new TH2D("dxvW2", "dx vs W2, all cuts; W2 (GeV); x_{HCAL}-x_{expect} (m)",binfac*W2fitmax, 0.0, W2fitmax, 600, dx_low,dx_high);
  TH1D *hdx_nocut = new TH1D( "dx_nocut", "HCal dx (m); m",hbinfac*hcal_fitrange ,hcalfit_low ,hcalfit_high);
  TH1D *hdx_cut = new TH1D( "dx_cut","HCal dx, all cuts; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );
  TH1D *hdx_cut_nofid = new TH1D( "dx_cut_nofid","HCal dx, all cuts but fiducial; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );
  TH1D *hdx_cut_failfid = new TH1D( "dx_cut_failfid","HCal dx, all cuts fail fiducial; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );
  TH1D *hdx_globcut = new TH1D( "dx_globcut","HCal dx, global cut; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );
  TH1D *hdx_glob_W2_cut = new TH1D( "dx_glob_W2_cut","HCal dx, global & W2 cuts; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );
  TH1D *hdx_nsigfid = new TH1D( "dx_nsigfid","HCal dx, nsigfid check; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );
  TH1D *hdx_nody = new TH1D( "hdx_cut_nody","HCal dx, all cuts but dy; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );

  TH1D *hdy_nocut = new TH1D( "dy_nocut", "HCal dy (m), no cuts; m", 250, dy_low, dy_high );
  TH1D *hdy_globcut = new TH1D( "dy_globcut", "HCal dy (m), global cut; m", 250, dy_low, dy_high );
  TH1D *hdy_glob_W2_cut = new TH1D( "dy_glob_W2_cut", "HCal dy (m), global & W2 cuts; m", 250, dy_low, dy_high );
  TH1D *hdy_cut_nofid = new TH1D( "dy_cut_nofid","HCal dy, all cuts but fiducial; y_{HCAL}-y_{expect} (m)", 250, dy_low, dy_high );
  TH1D *hdy_cut = new TH1D( "dy_cut","HCal dy, all cuts; y_{HCAL}-y_{expect} (m)", 250, dy_low, dy_high );
  TH1D *hdy_nody = new TH1D( "hdy_cut_nody","HCal dy, all cuts but dy; y_{HCAL}-y_{expect} (m)", 250, dy_low, dy_high );

  TH1D *hcoin_nocut = new TH1D( "hcoin_nocut", "HCal ADCt - BBCal ADCt, no cuts; ns", 200, -30, 30 );
  TH1D *hcoin_glob_W2_cut = new TH1D( "hcoin_glob_W2_cut", "HCal ADCt - BBCal ADCt, global & W2 cuts; ns", 200, -30, 30 );
  TH1D *hcoin_cut = new TH1D( "hcoin_cut", "HCal ADCt - BBCal ADCt, cuts; ns", 200, -30, 30 );
  TH1D *hcoin_pclus_glob_W2_cut = new TH1D( "hcoin_pclus_glob_W2_cut", "HCal ADCt - BBCal ADCt, pclus,global & W2 cuts; ns", 200, -30, 30 );
  TH1D *hcoin_pclus_cut = new TH1D( "hcoin_pclus_cut", "HCal ADCt - BBCal ADCt,pclus, cuts; ns", 200, -30, 30 );

  //general
  TH1D *hMott_cs = new TH1D( "hMott_cs", "Mott Cross Section, no cut; (GeV/c)^{-2}", 200, 0, 0.0002 );
  TH2D *hdx_xhcal = new TH2D("hdx_xhcal","dx vs xhcal; xhcal; x_{HCAL}-x_{expect} (m)",600,-3.0,3.0,hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high);
  TH2D *hdx_pN = new TH2D("hdx_pN", "dx vs nucleon momentum p_{N};p_{N} (GeV);x_{HCAL}-x_{expect} (m)",300,0.88*p_nuc_centr,1.12*p_nuc_centr,hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high);
  TH2D *hdx_prodef_pN = new TH2D("hdx_prodef_pN", "dx + proton deflection vs nucleon momentum p_{N};p_{N} (GeV);x_{HCAL}-x_{expect} + proton deflection (m)",300,0.88*p_nuc_centr,1.12*p_nuc_centr,hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high);
  TH1D *hdx_prodef = new TH1D( "dx_prodef","HCal dx + proton deflection; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );
  TH1D *hdx_defcut = new TH1D( "dx_defcut","HCal dx ; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );
  TH1D *hprodef = new TH1D( "hprodef","proton deflection; (m)", 300, 0.85*dx_pn, 1.3*dx_pn );

  //Added for cut stability studies
  TH1D* h_nsigx_fid = new TH1D("h_nsigx_fid", "nsigx_fid",200,-20,20);
  TH1D* h_nsigy_fid = new TH1D("h_nsigy_fid", "nsigy_fid",200,-20,20);

  //For anti-coin background
  TH1D *hcoin_coin_anticut = new TH1D( "hcoin_coin_anticut", "HCal ADCt - BBCal ADCt, coin_anticut; ns", 200, -100, 120 );
  TH1D *hdx_coin_anticut = new TH1D( "dx_coin_anticut","HCal dx, coin_anticut; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );
  TH1D *hdx_coin_anticut_wide = new TH1D( "dx_coin_anticut_wide","HCal dx, coin_anticut; x_{HCAL}-x_{expect} (m)", 400, -3.0, 2.0 );

  //Anti W2 background
  TH1D *h_W2_anticut_low  = new TH1D( "W2_anticut_low", "W2 (GeV^{2}) anticut lower region; GeV^{2}", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *h_W2_anticut_high  = new TH1D( "W2_anticut_high", "W2 (GeV^{2}) anticut higher region; GeV^{2}", binfac*4.5, 0.0, 4.5 );
  TH1D *h_W2_anticut_morehigh  = new TH1D( "W2_anticut_morehigh", "W2 (GeV^{2}) anticut more  higher region; GeV^{2}", binfac*4.5, 0.0, 4.5 );
  TH1D *h_W2_anticut_both  = new TH1D( "W2_anticut_both", "W2 (GeV^{2}) anticut both region; GeV^{2}", binfac*4.5, 0.0, 4.5 );

  TH1D *hdx_W2_anticut_low = new TH1D( "dx_W2_anticut_low","HCal dx, W2_anticut_low; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );
  TH1D *hdx_W2_anticut_low_wide = new TH1D( "dx_W2_anticut_low_wide","HCal dx, W2_anticut_low; x_{HCAL}-x_{expect} (m)", 400, -3.0, 2.0 );
  TH1D *hdx_W2_anticut_high = new TH1D( "dx_W2_anticut_high","HCal dx, W2_anticut_high; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );
  TH1D *hdx_W2_anticut_high_wide = new TH1D( "dx_W2_anticut_high_wide","HCal dx, W2_anticut_high; x_{HCAL}-x_{expect} (m)", 400, -3.0, 2.0 );
  TH1D *hdx_W2_anticut_morehigh = new TH1D( "dx_W2_anticut_morehigh","HCal dx, W2_anticut_morehigh; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );
  TH1D *hdx_W2_anticut_morehigh_wide = new TH1D( "dx_W2_anticut_morehigh_wide","HCal dx, W2_anticut_morehigh; x_{HCAL}-x_{expect} (m)", 400, -3.0, 2.0 );
  TH1D *hdx_W2_anticut_both = new TH1D( "dx_W2_anticut_both","HCal dx, W2_anticut_both; x_{HCAL}-x_{expect} (m)", hbinfac*hcal_fitrange, hcalfit_low, hcalfit_high );
  TH1D *hdx_W2_anticut_both_wide = new TH1D( "dx_W2_anticut_both_wide","HCal dx, W2_anticut_both; x_{HCAL}-x_{expect} (m)", 400, -3.0, 2.0 );

  //For physics extraction
  TH1D *h_Ebeam_corr_cut = new TH1D("Ebeam_corr_cut", "Beam Energy (GeV), all cuts; GeV", 300,EBeam-1.5,EBeam+1.5);
  TH1D *h_E_eprime_cut = new TH1D("E_eprime_cut", "Mag. Scattered Electron Momentum (GeV), all cuts; GeV", 300, 0.0, 5.0);
  TH1D *h_etheta_cut = new TH1D("etheta_cut","Theta for Scattered Electron from Reconstruct Track (radians), all cuts;",200, 0.0, 1.0);
  TH1D *h_pN_cut = new TH1D("h_pN_cut", "nucleon momentum p_{N}, all cuts;p_{N} (GeV)",300,0.8*p_nuc_centr,1.2*p_nuc_centr);
  TH1D *h_Q2_cut = new TH1D( "Q2_cut", "Q2 (GeV^{2}) all cuts; GeV^{2}", 300, 0.0, 6.0 );

  //allocate memory at each run
  TChain *C = nullptr;
  
  //create output tree
  TTree *Parse = new TTree("Parse","Analysis Tree");

  //uncut output tree variables
  double dx_out;
  double dy_out;
  double xexp_out;
  double yexp_out;
  double xhcal_out;
  double yhcal_out;
  double W2_out;
  double nu_out;
  double tau_out;
  double epsilon_out;
  double pcorr_out;
  double mott_out;
  double ehcal_out;
  double BBtot_e_out;
  double BBsh_e_out;
  double BBps_e_out;
  double hcal_atime_out;
  double BBsh_atime_out;
  double BBps_atime_out;
  double BBgem_nhits_out;
  double BBgem_ngoodhits_out;
  double BBgem_chi2ndf_out;
  double BBtr_x_out;
  double BBtr_y_out;
  double BBtr_p_out;
  double BBtr_vz_out;
  double BB_E_over_p_out;
  double BBtr_th_out;
  double BBtr_ph_out;
  double BBtr_r_x_out;
  double BBtr_r_y_out;
  double BBtr_r_th_out;
  double BBtr_r_ph_out;
  double coin_mean_out;
  double coin_sigma_out;
  double dyO_p_out;
  double dyO_n_out; 
  double dysig_p_out;
  double dysig_n_out;
  double dysig_n_fac_out;
  double dysig_p_fac_out;
  double dxO_p_out;
  double dxO_n_out;
  double dxsig_p_out;
  double dxsig_n_out;
  double dxsig_n_fac_out;
  double dxsig_p_fac_out;
  double nsigx_fid_out;
  double nsigy_fid_out;
  double W2low_out;
  double W2high_out;
  int num_hcal_clusid_out;
  int nclus_hcal_out;
  int hcal_clus_blk_out;
  int BBtr_n_out;
  int passGlobal_out;
  int HCalON_out;
  int passW2_out;
  int passCoin_out;
  int passFid_out;
  int run_out;
  int mag_out;
  int BBsh_nclus_out;
  int BBsh_nblk_out;
  int BBps_nclus_out;
  int BBps_nblk_out;

  double hcal_sh_atime_diff_out;
  double proton_deflection_out;
  double p_central_out;
  double rowblkHCAL_out;
  double colblkHCAL_out;
  double nblkHCAL_out;
  double ehcal_tree_out;

  //For physics extraction
  double Ebeam_corr_out;
  double E_eprime_out;
  double etheta_out;
  double p_N_out;
  double Q2_out;

  //setup new output tree branches
  Parse->Branch("dx", &dx_out, "dx/D");
  Parse->Branch("dy", &dy_out, "dy/D");
  Parse->Branch("xexp", &xexp_out, "xexp/D");
  Parse->Branch("yexp", &yexp_out, "yexp/D");
  Parse->Branch("xhcal", &xhcal_out, "xhcal/D");
  Parse->Branch("yhcal", &yhcal_out, "yhcal/D");
  Parse->Branch("W2", &W2_out, "W2/D");
  Parse->Branch("nu", &nu_out, "nu/D");
  Parse->Branch("tau", &tau_out, "tau/D");
  Parse->Branch("epsilon", &epsilon_out, "epsilon/D");
  Parse->Branch("pcorr", &pcorr_out, "pcorr/D");
  Parse->Branch("mott", &mott_out, "mott/D");
  Parse->Branch("ehcal", &ehcal_out, "ehcal/D");
  Parse->Branch("ehcal_tree", &ehcal_tree_out, "ehcal_tree/D");
  Parse->Branch("BBtot_e", &BBtot_e_out, "BBtot_e/D");
  Parse->Branch("BBsh_e", &BBsh_e_out, "BBsh_e/D");
  Parse->Branch("BBsh_nclus", &BBsh_nclus_out, "BBsh_nclus/I");
  Parse->Branch("BBsh_nblk", &BBsh_nblk_out, "BBsh_nblk/I");
  Parse->Branch("BBps_e", &BBps_e_out, "BBps_e/D");
  Parse->Branch("BBps_nclus", &BBps_nclus_out, "BBps_nclus/I");
  Parse->Branch("BBps_nblk", &BBps_nblk_out, "BBps_nblk/I");
  Parse->Branch("hcal_atime", &hcal_atime_out, "hcal_atime/D");
  Parse->Branch("BBsh_atime", &BBsh_atime_out, "BBsh_atime/D");
  Parse->Branch("BBps_atime", &BBps_atime_out, "BBps_atime/D");
  Parse->Branch("BBgem_nhits", &BBgem_nhits_out, "BBgem_nhits/D");
  Parse->Branch("BBgem_ngoodhits", &BBgem_ngoodhits_out, "BBgem_ngoodhits/D");
  Parse->Branch("BBgem_chi2ndf", &BBgem_chi2ndf_out, "BBgem_chi2ndf/D");
  Parse->Branch("BBtr_x", &BBtr_x_out, "BBtr_x/D");
  Parse->Branch("BBtr_y", &BBtr_y_out, "BBtr_y/D");
  Parse->Branch("BBtr_p", &BBtr_p_out, "BBtr_p/D");
  Parse->Branch("BBtr_vz", &BBtr_vz_out, "BBtr_vz/D");
  Parse->Branch("BB_E_over_p", &BB_E_over_p_out, "BB_E_over_p/D");
  Parse->Branch("BBtr_th", &BBtr_th_out, "BBtr_th/D");
  Parse->Branch("BBtr_ph", &BBtr_ph_out, "BBtr_ph/D");
  Parse->Branch("BBtr_r_x", &BBtr_r_x_out, "BBtr_r_x/D");
  Parse->Branch("BBtr_r_y", &BBtr_r_y_out, "BBtr_r_y/D");
  Parse->Branch("BBtr_r_th", &BBtr_r_th_out, "BBtr_r_th/D");
  Parse->Branch("BBtr_r_ph", &BBtr_r_ph_out, "BBtr_r_ph/D");
  Parse->Branch("coin_mean", &coin_mean_out, "coin_mean/D");
  Parse->Branch("coin_sigma", &coin_sigma_out, "coin_sigma/D");
  Parse->Branch("dyO_p", &dyO_p_out, "dyO_p/D");
  Parse->Branch("dyO_n", &dyO_n_out, "dyO_n/D");
  Parse->Branch("dysig_p", &dysig_p_out, "dysig_p/D");
  Parse->Branch("dysig_n", &dysig_n_out, "dysig_n/D");
  Parse->Branch("dysig_n_fac", &dysig_n_fac_out, "dysig_n_fac/D");
  Parse->Branch("dysig_p_fac", &dysig_p_fac_out, "dysig_p_fac/D");
  Parse->Branch("dxO_p", &dxO_p_out, "dxO_p/D");
  Parse->Branch("dxO_n", &dxO_n_out, "dxO_n/D");
  Parse->Branch("dxsig_p", &dxsig_p_out, "dxsig_p/D");
  Parse->Branch("dxsig_n", &dxsig_n_out, "dxsig_n/D");
  Parse->Branch("dxsig_n_fac", &dxsig_n_fac_out, "dxsig_n_fac/D");
  Parse->Branch("dxsig_p_fac", &dxsig_p_fac_out, "dxsig_p_fac/D");
  Parse->Branch("nsigx_fid", &nsigx_fid_out , "nsigx_fid/D");
  Parse->Branch("nsigy_fid", &nsigy_fid_out , "nsigy_fid/D");
  Parse->Branch("W2low", &W2low_out, "W2low/D");
  Parse->Branch("W2high", &W2high_out, "W2high/D");
  Parse->Branch("hcal_sh_atime_diff", &hcal_sh_atime_diff_out, "hcal_sh_atime_diff/D");
  Parse->Branch("proton_deflection", &proton_deflection_out, "proton_deflection/D");
  Parse->Branch("p_central", &p_central_out, "p_centrial/D");

  Parse->Branch("num_hcal_clusid", &num_hcal_clusid_out, "num_hcal_clusid/I");
  Parse->Branch("hcal_clus_blk", &hcal_clus_blk_out, "hcal_clus_blk/I");
  Parse->Branch("nclus_hcal", &nclus_hcal_out,"nclus_hcal/I");
  Parse->Branch("BBtr_n", &BBtr_n_out, "BBtr_n/I");
  Parse->Branch("passGlobal", &passGlobal_out,"passGlobal/I");
  Parse->Branch("HCalON", &HCalON_out,"HCalON/I");
  Parse->Branch("passW2", &passW2_out,"passW2/I");
  Parse->Branch("passCoin", &passCoin_out,"passCoin/I");
  Parse->Branch("passFid", &passFid_out, "passFid/I");
  Parse->Branch("run", &run_out, "run/I");
  Parse->Branch("mag", &mag_out, "mag/I");

  Parse->Branch( "nblkHCAL", &nblkHCAL_out, "nblkHCAL/D" );
  Parse->Branch( "rowblkHCAL",&rowblkHCAL_out, "rowblkHCAL/D" );
  Parse->Branch( "colblkHCAL", &colblkHCAL_out, "colblkHCAL/D" );

  //For physics extraction
  Parse->Branch("Ebeam_corr", &Ebeam_corr_out, "Ebeam_corr/D");
  Parse->Branch("E_eprime", &E_eprime_out, "E_eprime/D");
  Parse->Branch("etheta", &etheta_out, "etheta/D");
  Parse->Branch("p_N", &p_N_out, "p_N/D");
  Parse->Branch("Q2", &Q2_out, "Q2/D");

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
  double x_hcal,y_hcal,e_hcal,nclus_hcal,idx_hcal, nblkHCAL;
  
  C->SetBranchStatus("sbs.hcal.nblk",1);
  C->SetBranchStatus("sbs.hcal.x",1);
  C->SetBranchStatus("sbs.hcal.y",1);
  C->SetBranchStatus("sbs.hcal.e",1);
  C->SetBranchStatus("sbs.hcal.nclus",1); 
  C->SetBranchStatus("sbs.hcal.index",1);

  C->SetBranchAddress("sbs.hcal.nblk",&nblkHCAL);
  C->SetBranchAddress("sbs.hcal.x", &x_hcal);
  C->SetBranchAddress("sbs.hcal.y", &y_hcal);
  C->SetBranchAddress("sbs.hcal.e", &e_hcal);
  C->SetBranchAddress("sbs.hcal.nclus", &nclus_hcal);
  C->SetBranchAddress("sbs.hcal.index", &idx_hcal);

  //HCal cluster branches
  double hcal_clus_atime[exp_constants::maxclus], hcal_clus_e[exp_constants::maxclus], hcal_clus_x[exp_constants::maxclus], hcal_clus_y[exp_constants::maxclus],hcal_clus_nblk[exp_constants::maxclus], rowblkHCAL[exp_constants::maxclus], colblkHCAL[exp_constants::maxclus];
  int num_hcal_clusid;

  C->SetBranchStatus("sbs.hcal.clus.atime",1);
  C->SetBranchStatus("sbs.hcal.clus.e",1);
  C->SetBranchStatus("sbs.hcal.clus.x",1);
  C->SetBranchStatus("sbs.hcal.clus.y",1);
  C->SetBranchStatus("sbs.hcal.clus.e",1);
  C->SetBranchStatus("Ndata.sbs.hcal.clus.id", 1);
  C->SetBranchStatus("sbs.hcal.clus.nblk", 1);
  C->SetBranchStatus("sbs.hcal.clus_blk.row", 1);
  C->SetBranchStatus("sbs.hcal.clus_blk.col", 1);

  C->SetBranchAddress("sbs.hcal.clus.atime", &hcal_clus_atime);
  C->SetBranchAddress("sbs.hcal.clus.e", &hcal_clus_e);
  C->SetBranchAddress("sbs.hcal.clus.x", &hcal_clus_x);
  C->SetBranchAddress("sbs.hcal.clus.y", &hcal_clus_y);
  C->SetBranchAddress("Ndata.sbs.hcal.clus.id", &num_hcal_clusid);
  C->SetBranchAddress("sbs.hcal.clus.nblk", &hcal_clus_nblk);

  C->SetBranchAddress("sbs.hcal.clus_blk.row",&rowblkHCAL);
  C->SetBranchAddress("sbs.hcal.clus_blk.col",&colblkHCAL);

  //BBCal shower
  double atime_sh, e_sh, nclus_sh, nblk_sh; 

  C->SetBranchStatus( "bb.sh.atimeblk", 1 );
  C->SetBranchStatus( "bb.sh.e", 1 );
  C->SetBranchStatus( "bb.sh.nclus", 1 );
  C->SetBranchStatus( "bb.sh.nblk", 1 );

  C->SetBranchAddress("bb.sh.atimeblk", &atime_sh);
  C->SetBranchAddress("bb.sh.e", &e_sh);  
  C->SetBranchAddress("bb.sh.nclus", &nclus_sh);
  C->SetBranchAddress("bb.sh.nblk", &nblk_sh);


  //BBCal preshower
  
  double atime_ps, e_ps, nclus_ps, nblk_ps;

  C->SetBranchStatus( "bb.ps.atimeblk", 1 );
  C->SetBranchStatus( "bb.ps.e", 1 );
  C->SetBranchStatus( "bb.ps.nclus", 1 );
  C->SetBranchStatus( "bb.ps.nblk", 1 );

  C->SetBranchAddress("bb.ps.atimeblk", &atime_ps);       
  C->SetBranchAddress("bb.ps.e", &e_ps);
  C->SetBranchAddress("bb.ps.nclus", &nclus_ps);
  C->SetBranchAddress("bb.ps.nblk", &nblk_ps);

  //BBGEM hits

  double gem_hits[maxtracks],gem_goodhits[maxtracks],gem_ChiSqr[maxtracks];  

  C->SetBranchStatus("bb.gem.track.nhits", 1);
  C->SetBranchStatus("bb.gem.track.ngoodhits",1);
  C->SetBranchStatus("bb.gem.track.chi2ndf",1);

  C->SetBranchAddress("bb.gem.track.nhits",&gem_hits);
  C->SetBranchAddress("bb.gem.track.ngoodhits",&gem_goodhits);
  C->SetBranchAddress("bb.gem.track.chi2ndf",&gem_ChiSqr);
  // track branches

  double ntrack, tr_px[maxtracks], tr_py[maxtracks], tr_pz[maxtracks], tr_p[maxtracks], tr_x[maxtracks], tr_y[maxtracks], tr_vx[maxtracks], tr_vy[maxtracks], tr_vz[maxtracks], tr_r_x[maxtracks], tr_r_y[maxtracks],tr_r_th[maxtracks], tr_r_ph[maxtracks], tr_th[maxtracks], tr_ph[maxtracks];


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

  C->SetBranchAddress("bb.tr.n",&ntrack);
  C->SetBranchAddress("bb.tr.px",&tr_px);
  C->SetBranchAddress("bb.tr.py",&tr_py);
  C->SetBranchAddress("bb.tr.pz",&tr_pz);
  C->SetBranchAddress("bb.tr.p",&tr_p);
  C->SetBranchAddress("bb.tr.x",&tr_x);
  C->SetBranchAddress("bb.tr.y",&tr_y);
  C->SetBranchAddress("bb.tr.vx",&tr_vx);
  C->SetBranchAddress("bb.tr.vy",&tr_vy);
  C->SetBranchAddress("bb.tr.vz",&tr_vz);
  C->SetBranchAddress("bb.tr.th",&tr_th);
  C->SetBranchAddress("bb.tr.ph",&tr_ph);
  C->SetBranchAddress("bb.tr.r_th",&tr_r_th);
  C->SetBranchAddress("bb.tr.r_ph",&tr_r_ph);
  C->SetBranchAddress("bb.tr.r_x",&tr_r_x);
  C->SetBranchAddress("bb.tr.r_y",&tr_r_y);

  //ekine branches

  double ekine_Q2, ekine_W2, ekine_eps, ekine_nu, ekine_qx, ekine_qy, ekine_qz;

  C->SetBranchStatus("e.kine.Q2",1);
  C->SetBranchStatus("e.kine.W2",1);
  C->SetBranchStatus("e.kine.epsilon",1);
  C->SetBranchStatus("e.kine.nu",1);
  C->SetBranchStatus("e.kine.q_x",1);
  C->SetBranchStatus("e.kine.q_y",1);
  C->SetBranchStatus("e.kine.q_z",1);

  C->SetBranchAddress("e.kine.Q2", &ekine_Q2);
  C->SetBranchAddress("e.kine.W2", &ekine_W2); 
  C->SetBranchAddress("e.kine.epsilon", &ekine_eps);
  C->SetBranchAddress("e.kine.nu", &ekine_nu);
  C->SetBranchAddress("e.kine.q_x", &ekine_qx);
  C->SetBranchAddress("e.kine.q_y", &ekine_qy);
  C->SetBranchAddress("e.kine.q_z", &ekine_qz);

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
    	}
        //Is true if failed global cut
    	bool failglobal = cuts::failedGlobal(GlobalCut);
	
	///////////
	//Electron-arm physics calculations
	
	//BB E/p calculated. Should be the same as the tree variable. But it does not cause seg faults
	double BB_E_over_p = (e_sh+e_ps)/tr_p[0];

	//correct beam energy from vertex information
	double Eloss = physics::getEloss(tr_vz[0],target);
	double Ecorr = physics::getEcorr(Ebeam,Eloss);

	//make the vertex, assuming only beam direction
	TVector3 vertex = physics::getVertex(tr_vz[0]);

	//reconstructed momentum from track momentum information, corrected for mean energy loss. Still need to include losses from Al shielding or target windows later
	double pcorr = physics::getp_recon_corr(tr_p[0],Eloss_outgoing); 

        //four momentum vector for electron beam with correted Energy value
        TLorentzVector pbeam = physics::getpBeam(Ecorr);

	//four momentum for scattered electron based on reconstruction
	TLorentzVector p_eprime = physics::getp_eprime(tr_px[0],tr_py[0],tr_pz[0],tr_p[0],pcorr);

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
		Q2 = physics::getQ2(ekine_Q2);
		W2 = physics::getW2(ekine_W2);
		nu = physics::getnu(ekine_nu);
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
	int intime_idx = physics::cluster_intime_select(num_hcal_clusid,hcal_clus_atime,atime_sh,hcal_clus_e,coin_mean,coin_sig_fac,coin_profile_sig,hcalemin);
	
	//Assume that the itime analysis is sufficient to find the best cluster in HCal
	int clus_idx_best = intime_idx;

	//calculate important information from best cluster
	double xhcal_bestclus = hcal_clus_x[clus_idx_best];
 	double yhcal_bestclus = hcal_clus_y[clus_idx_best];
	double xhcal_pclus = hcal_clus_x[0];
        double yhcal_pclus = hcal_clus_y[0];
	double dx_bestclus = physics::get_dx(xhcal_bestclus,xhcal_expect);
	double dx_pclus = physics::get_dx(xhcal_pclus,xhcal_expect);
 	double dy_bestclus = physics::get_dy(yhcal_bestclus,yhcal_expect);	
	double hcal_atime_bestclus = hcal_clus_atime[clus_idx_best];
	double coin_bestclus = hcal_atime_bestclus - atime_sh;
	double coin_pclus = hcal_clus_atime[0] - atime_sh;
	double hcal_e_bestclus = hcal_clus_e[clus_idx_best];
	int hcal_nblk_bestclus = (int) hcal_clus_nblk[clus_idx_best];	

	//calculate the number of sigma away from fiducial boundaries. Store info for later
        double nsigx_fid = cuts::calculate_nsigma_fid_x(xhcal_expect,dxsig_p,dxsig_n,dx_pn,hcalaa);
        double nsigy_fid = cuts::calculate_nsigma_fid_y(yhcal_expect,dysig_p,hcalaa);

	
	//setup booleans for cuts later. Save boolean values to tree
	//global is above
	
	//HCal active area
	bool hcalaa_ON = cuts::hcalaa_ON(xhcal_bestclus,yhcal_bestclus,hcalaa);
	bool hcalaa_ON_exp = cuts::hcalaa_ON(xhcal_expect,yhcal_expect,hcalaa);

	//W2 elastic boolean
	bool goodW2 = cuts::goodW2(W2,W2_low,W2_high);
	bool antiW2_low = W2 > 0.0 && W2 < W2_low;
	bool antiW2_high = W2 > W2_high;
	bool antiW2_morehigh = W2 > (W2_high+1.1);

	//good dy boolean
	bool good_dy = cuts::good_dy(dy_bestclus,dyO_p,dysig_cut_fac,dysig_cut);

	//good coincidence time cut
	bool passCoin = cuts::passCoin(coin_bestclus,coin_mean,coin_sig_fac,coin_sig);
	//Used for coin anti cut
	bool passCoin_pclus = cuts::passCoin(coin_pclus,coin_mean,coin_sig_fac,coin_sig);
	
	//good fiducial cut
	bool passFid = cuts::hcalfid_IN(xhcal_expect,yhcal_expect,dx_pn,hcalfid);	

	//pass HCal E
	bool passHCalE = cuts::passHCalE(hcal_e_bestclus,hcalemin);

	//pass HCal num clus
	bool passHCal_Nclus = cuts::passHCal_NClus(nclus_hcal,hcalnclusmin);

	//pass NSig Fid check
	bool passNSigFid = cuts::passNsigFid(nsigx_fid,nsigy_fid);

	//Fill analysis tree variables before making cuts
	dx_out = dx_bestclus;
 	dy_out = dy_bestclus;
  	xexp_out = xhcal_expect;
  	yexp_out = yhcal_expect;
  	xhcal_out = xhcal_bestclus;
 	yhcal_out = yhcal_bestclus;
  	W2_out = W2;
  	nu_out = nu;
	tau_out = tau;
	epsilon_out = epsilon;
	pcorr_out = pcorr;
	mott_out = Mott_CS;
  	ehcal_out = hcal_e_bestclus;
  	ehcal_tree_out = e_hcal;
	BBtot_e_out = e_sh+e_ps;
  	BBsh_e_out = e_sh;
  	BBps_e_out = e_ps;
  	hcal_atime_out = hcal_atime_bestclus;
  	BBsh_atime_out = atime_sh;
  	BBps_atime_out = atime_ps;
  	BBgem_nhits_out = gem_hits[0];
  	BBgem_ngoodhits_out = gem_goodhits[0];
	BBgem_chi2ndf_out = gem_ChiSqr[0];
	BBtr_x_out = tr_x[0];
  	BBtr_y_out = tr_y[0];
  	BBtr_p_out = tr_p[0];
  	BBtr_vz_out = tr_vz[0];
	BB_E_over_p_out = BB_E_over_p;
	BBtr_th_out = tr_th[0];
	BBtr_ph_out = tr_ph[0];
	BBtr_r_x_out = tr_r_x[0];
	BBtr_r_y_out = tr_r_y[0];
	BBtr_r_th_out = tr_r_th[0];
	BBtr_r_ph_out = tr_r_ph[0];
	coin_mean_out = coin_mean;
  	coin_sigma_out = coin_sig;
  	dyO_p_out = dyO_p;
  	dyO_n_out = dyO_n;
  	dysig_p_out = dysig_p;
  	dysig_n_out = dysig_n;
  	dysig_n_fac_out = dysig_n_fac;
  	dysig_p_fac_out = dysig_p_fac;
  	dxO_p_out = dxO_p;
  	dxO_n_out = dxO_n;
  	dxsig_p_out = dxsig_p;
  	dxsig_n_out = dxsig_n;
  	dxsig_n_fac_out = dxsig_n_fac;
 	dxsig_p_fac_out = dxsig_p_fac;
  	nsigx_fid_out = nsigx_fid;
	nsigy_fid_out = nsigy_fid;
	W2low_out = W2_low;	
	W2high_out = W2_high;
	num_hcal_clusid_out = num_hcal_clusid ;
 	hcal_clus_blk_out = hcal_nblk_bestclus;
  	nclus_hcal_out = nclus_hcal;
	BBsh_nclus_out = (int) nclus_sh;
	BBsh_nblk_out = (int) nblk_sh;
	BBps_nclus_out = (int) nclus_ps;
        BBps_nblk_out = (int) nblk_ps;
	BBtr_n_out = ntrack;
  	passGlobal_out = (int) !failglobal;
	HCalON_out = (int) hcalaa_ON;
	passW2_out = (int) goodW2;
	passCoin_out = (int) passCoin;
	passFid_out = (int) passFid;
  	run_out = run ;
  	mag_out = field;
	hcal_sh_atime_diff_out = coin_bestclus;
	proton_deflection_out = proton_deflection;
	p_central_out = pcentral;
	nblkHCAL_out = nblkHCAL;
	rowblkHCAL_out = rowblkHCAL[clus_idx_best];
	colblkHCAL_out = colblkHCAL[clus_idx_best];

	//For physics extraction
	Ebeam_corr_out = Ecorr;
	E_eprime_out =  p_eprime.E();
	etheta_out = etheta;
	p_N_out = p_N.Vect().Mag();
	Q2_out = Q2;

	//Fill histograms of global cut parameters here without any restrictions
	h_ntracks->Fill(ntrack);	
	h_PS_E->Fill(e_ps);
	h_vert_z->Fill(tr_vz[0]);
	h_HCal_E->Fill(e_hcal);
        h_HCal_E_best->Fill(hcal_e_bestclus);
	h_HCal_nclus->Fill(nclus_hcal);
	h_nhits->Fill(gem_hits[0]);
	h_bbtrp_nocut->Fill(tr_p[0]);
	h_SH_nclus->Fill(nclus_sh);
	h_bbEoverp_nocut->Fill(BB_E_over_p);
	h_optics_xdir_nocut->Fill(tr_r_x[0]-0.9*tr_r_th[0]);
	h_W2_optics_xdir_nocut->Fill(tr_r_x[0]-0.9*tr_r_th[0],W2);
	h_optics_ydir_nocut->Fill(tr_r_y[0]-0.9*tr_r_ph[0]);
        h_W2_optics_ydir_nocut->Fill(tr_r_y[0]-0.9*tr_r_ph[0],W2);
	h_track_chi2ndf_nocut->Fill(gem_ChiSqr[0]);
	h_xexp_optics_xdir_nocut->Fill(tr_r_x[0]-0.9*tr_r_th[0],xhcal_expect);
	h_yexp_optics_ydir_nocut->Fill(tr_r_y[0]-0.9*tr_r_ph[0],yhcal_expect);

	hdx_nocut->Fill(dx_bestclus);
	hcoin_nocut->Fill(coin_bestclus);
	hdy_nocut->Fill(dy_bestclus);
	hxy_expect_nocut->Fill(yhcal_expect,xhcal_expect);
	hx_expect_nocut->Fill(xhcal_expect);
	hy_expect_nocut->Fill(yhcal_expect);
	h_W2_nocut->Fill(W2);

	hxy_expect_nocut_n->Fill(yhcal_expect,xhcal_expect);
	hxy_expect_nocut_p->Fill(yhcal_expect,xhcal_expect-dx_pn);

	//Fill some histograms here after basic global cuts
	if(!failglobal){
	//global parameter checks
	h_ntracks_globcut->Fill(ntrack);
        h_PS_E_globcut->Fill(e_ps);
        h_vert_z_globcut->Fill(tr_vz[0]);
        h_HCal_E_globcut->Fill(e_hcal);
        h_HCal_E_best_globcut->Fill(hcal_e_bestclus);
	h_HCal_nclus_globcut->Fill(nclus_hcal);
        h_nhits_globcut->Fill(gem_hits[0]);
        h_bbtrp_globcut->Fill(tr_p[0]);
        h_SH_nclus_globcut->Fill(nclus_sh);
	h_TPS_SH_globcut->Fill(e_ps+e_sh);
	h_bbEoverp_globcut->Fill(BB_E_over_p);
	h_optics_xdir_globcut->Fill(tr_r_x[0]-0.9*tr_r_th[0]);
        h_W2_optics_xdir_globcut->Fill(tr_r_x[0]-0.9*tr_r_th[0],W2);
        h_optics_ydir_globcut->Fill(tr_r_y[0]-0.9*tr_r_ph[0]);
        h_W2_optics_ydir_globcut->Fill(tr_r_y[0]-0.9*tr_r_ph[0],W2);
	h_track_chi2ndf_globcut->Fill(gem_ChiSqr[0]);
	h_xexp_optics_xdir_globcut->Fill(tr_r_x[0]-0.9*tr_r_th[0],xhcal_expect);
        h_yexp_optics_ydir_globcut->Fill(tr_r_y[0]-0.9*tr_r_ph[0],yhcal_expect);

	//physics quantities
	hxy_globcut->Fill(yhcal_bestclus,xhcal_bestclus);
	hx_expect_globcut->Fill(xhcal_expect);
        hy_expect_globcut->Fill(yhcal_expect);
	h_W2_globcut->Fill(W2);
	h_Q2_globcut->Fill(Q2);
	hxy_expect_globcut->Fill(yhcal_expect,xhcal_expect);
	hdxdy_globcut->Fill(dy_bestclus, dx_bestclus);
	hdx_globcut->Fill(dx_bestclus);
	hdy_globcut->Fill(dy_bestclus);
	
	}
	
	//Fill some histograms if pass global cut and W2 cut. Mostly just e-arm cuts
	if(!failglobal && goodW2){
	hxy_glob_W2_cut->Fill(yhcal_bestclus,xhcal_bestclus);
        h_W2_glob_W2_cut->Fill(W2);
        hxy_expect_glob_W2_cut->Fill(yhcal_expect,xhcal_expect);
        hdxdy_glob_W2_cut->Fill(dy_bestclus, dx_bestclus);
        hdx_glob_W2_cut->Fill(dx_bestclus);
        hdy_glob_W2_cut->Fill(dy_bestclus);
	hcoin_glob_W2_cut->Fill(coin_bestclus);
	hcoin_pclus_glob_W2_cut->Fill(coin_pclus);
	hxy_expect_n->Fill(yhcal_expect,xhcal_expect);
	hxy_expect_p->Fill(yhcal_expect,(xhcal_expect-dx_pn));
	hMott_cs->Fill(Mott_CS);

	hxy_expect_glob_W2_cut_n->Fill(yhcal_expect,xhcal_expect);
	hxy_expect_glob_W2_cut_p->Fill(yhcal_expect,xhcal_expect-dx_pn);
	}
	
	//Now let's add in our first major hadron arm cut along with all the cuts from before.
	if(!failglobal && goodW2 && hcalaa_ON){
	hxy_acceptancecut->Fill(yhcal_bestclus,xhcal_bestclus);
	}

	//related checks for proton deflection
	if(!failglobal && passHCalE && passHCal_Nclus && goodW2 && good_dy){
	hdx_xhcal->Fill(xhcal_bestclus,dx_bestclus);
	hdx_pN->Fill(p_N.Vect().Mag(),dx_bestclus);
	hdx_prodef_pN->Fill(p_N.Vect().Mag(),dx_bestclus+proton_deflection);
	hdx_prodef->Fill(dx_bestclus+proton_deflection);
	hdx_defcut->Fill(dx_bestclus);
	hprodef->Fill(proton_deflection);
	}

	//e-arm and h-arm cuts, except fiducial
	if(!failglobal && passHCalE && passHCal_Nclus && goodW2 && hcalaa_ON && passCoin && good_dy ){
	hdxdy_nofid->Fill(dy_bestclus, dx_bestclus);
	hdx_cut_nofid->Fill(dx_bestclus);
	hdy_cut_nofid->Fill(dy_bestclus);
	hdx_xhcal->Fill(xhcal_bestclus,dx_bestclus);
		if(!passFid){
		hdx_cut_failfid->Fill(dx_bestclus);
		}

	}	
	
	//all cuts but W2
	if(!failglobal && passHCalE && passHCal_Nclus && hcalaa_ON && passCoin && good_dy && passFid){
	h_W2_notW2_cut->Fill(W2);

	}
	
	//Essentially all cuts and anti coin cut
	if(!failglobal && passHCalE && passHCal_Nclus && goodW2 && hcalaa_ON && !passCoin_pclus && good_dy && passFid){
	hcoin_coin_anticut->Fill(coin_pclus);
	hdx_coin_anticut->Fill(dx_pclus);
	hdx_coin_anticut_wide->Fill(dx_pclus);
	}

	//Essentially all cuts and anti W2 cut low region
	if(!failglobal && passHCalE && passHCal_Nclus && antiW2_low && hcalaa_ON && passCoin && good_dy && passFid){
	h_W2_anticut_low->Fill(W2);
	hdx_W2_anticut_low->Fill(dx_bestclus);
	hdx_W2_anticut_low_wide->Fill(dx_bestclus);
	}
	//Essentially all cuts and anti W2 cut high region
	if(!failglobal && passHCalE && passHCal_Nclus && antiW2_high && hcalaa_ON && passCoin && good_dy && passFid){
        h_W2_anticut_high->Fill(W2);
        hdx_W2_anticut_high->Fill(dx_bestclus);
        hdx_W2_anticut_high_wide->Fill(dx_bestclus);
        }
	
	//Essentially all cuts and anti W2 cut high region
        if(!failglobal && passHCalE && passHCal_Nclus && antiW2_morehigh && hcalaa_ON && passCoin && good_dy && passFid){
        h_W2_anticut_morehigh->Fill(W2);
        hdx_W2_anticut_morehigh->Fill(dx_bestclus);
        hdx_W2_anticut_morehigh_wide->Fill(dx_bestclus);
        }

	//Essentially all cuts and anti W2 cut both region
	if(!failglobal && passHCalE && passHCal_Nclus && (antiW2_low || antiW2_high) && hcalaa_ON && passCoin && good_dy && passFid){
        h_W2_anticut_both->Fill(W2);
        hdx_W2_anticut_both->Fill(dx_bestclus);
        hdx_W2_anticut_both_wide->Fill(dx_bestclus);
        }
	
	//all cuts but HCal E
	if(!failglobal  && passHCal_Nclus && goodW2 && hcalaa_ON && passCoin && good_dy && passFid){
	h_HCal_E_best_cut_noHCalE->Fill(hcal_e_bestclus);
	h_HCal_E_cut_noHCalE->Fill(e_hcal);
	
	}
	
	
	//Essentially all cuts but dy (visualization)
	if(!failglobal && passHCalE && passHCal_Nclus && goodW2 && hcalaa_ON && passCoin  && passFid){
	hdxdy_nody->Fill(dy_bestclus, dx_bestclus);
	hdx_nody->Fill(dx_bestclus);
	hdy_nody->Fill(dy_bestclus);
	}
	
	//all cuts
	if(!failglobal && passHCalE && passHCal_Nclus && goodW2 && hcalaa_ON && passCoin && good_dy && passFid){
	h_ntracks_cut->Fill(ntrack);
        h_PS_E_cut->Fill(e_ps);
        h_vert_z_cut->Fill(tr_vz[0]);
        h_HCal_E_cut->Fill(e_hcal);
        h_HCal_E_best_cut->Fill(hcal_e_bestclus);
	h_HCal_nclus_cut->Fill(nclus_hcal);
        h_nhits_cut->Fill(gem_hits[0]);
        h_bbtrp_cut->Fill(tr_p[0]);
        h_SH_nclus_cut->Fill(nclus_sh);
        h_TPS_SH_cut->Fill(e_ps+e_sh);
        h_bbEoverp_cut->Fill(BB_E_over_p);
	hxy_cut->Fill(yhcal_bestclus,xhcal_bestclus);
        h_W2_cut->Fill(W2);
        hxy_expect_cut->Fill(yhcal_expect,xhcal_expect);
        hdxdy_cut->Fill(dy_bestclus, dx_bestclus);
        hdx_cut->Fill(dx_bestclus);
        hdy_cut->Fill(dy_bestclus);
	hcoin_cut->Fill(coin_bestclus);
        hcoin_pclus_cut->Fill(coin_pclus);
        hxy_expect_fidcutn->Fill(yhcal_expect,xhcal_expect);
	h_optics_xdir_cut->Fill(tr_r_x[0]-0.9*tr_r_th[0]);
        h_W2_optics_xdir_cut->Fill(tr_r_x[0]-0.9*tr_r_th[0],W2);
        h_optics_ydir_cut->Fill(tr_r_y[0]-0.9*tr_r_ph[0]);
        h_W2_optics_ydir_cut->Fill(tr_r_y[0]-0.9*tr_r_ph[0],W2);
	h_track_chi2ndf_cut->Fill(gem_ChiSqr[0]);
	h_xexp_optics_xdir_cut->Fill(tr_r_x[0]-0.9*tr_r_th[0],xhcal_expect);
        h_yexp_optics_ydir_cut->Fill(tr_r_y[0]-0.9*tr_r_ph[0],yhcal_expect);
	hx_expect_cut->Fill(xhcal_expect);
        hy_expect_cut->Fill(yhcal_expect);

	hxy_expect_fidcutp->Fill(yhcal_expect,(xhcal_expect-dx_pn));
	hdxvE->Fill(hcal_e_bestclus,dx_bestclus);
	hdxvW2->Fill(W2,dx_bestclus);
	
	//For physics extraction
	h_Ebeam_corr_cut->Fill(Ecorr);
	h_E_eprime_cut->Fill(p_eprime.E());
	h_etheta_cut->Fill(etheta);
	h_pN_cut->Fill(p_N.Vect().Mag());
	h_Q2_cut->Fill(Q2);

	}
	
	if(!failglobal && passHCalE && passHCal_Nclus && goodW2 && hcalaa_ON && passCoin && good_dy && !passFid){
	hxy_expect_failedfid->Fill(yhcal_expect,xhcal_expect);
	}

	//For cut stability
	h_nsigx_fid ->Fill(nsigx_fid);
	h_nsigy_fid ->Fill(nsigy_fid);
	if(passNSigFid){
	hdx_nsigfid->Fill(dx_bestclus);
	}

	
	//Fill the analysis tree
	Parse->Fill();
	
	}//end event loop

	// reset chain for the next run config
	C->Reset();
  
  }//end loop over the run numbers
  
  //make lines for active area on HCal
  vector<TLine*> Lines_aa = plots::setupLines(hcalaa,2,kRed);
 
 //make lines for fiducial region
  vector<TLine*> Lines_Fid = plots::setupLines(hcalfid,2,kMagenta);
  TLine *LineFidPro = plots::setupLine_Horz(hcalfid[2],hcalfid[3],hcalfid[0]+dx_pn,2,kMagenta,2);

  //make lines for physical HCal position
  vector<TLine*> Lines_pos = plots::setupLines(hcalpos,2,kGreen);

  TCanvas* c0 = plots::plotAcceptance_Check("c0",Lines_pos,Lines_aa,Lines_Fid,hxy_globcut,hxy_acceptancecut);
  TCanvas* c1 = plots::plotFid_Check("c1",Lines_pos,Lines_aa,Lines_Fid,LineFidPro,hxy_expect_glob_W2_cut,hxy_expect_cut,hxy_expect_failedfid);
  TCanvas* c2 = plots::plotFid_Hypothesis_Check("c2",Lines_pos,Lines_aa,Lines_Fid,LineFidPro,hxy_expect_fidcutn,hxy_expect_fidcutp);
  TCanvas* c3 = plots::plotFid_Hypothesis_Check("c3",Lines_pos,Lines_aa,Lines_Fid,LineFidPro,hxy_expect_nocut_n,hxy_expect_nocut_p);
  TCanvas* c4 = plots::plotFid_Hypothesis_Check("c4",Lines_pos,Lines_aa,Lines_Fid,LineFidPro,hxy_expect_glob_W2_cut_n,hxy_expect_glob_W2_cut_p);

  //Write stuff to a pdf
  TString plotname = outfile;
  plotname.ReplaceAll(".root",".pdf");
  TString start = Form("%s%s",plotname.Data(),"(");
  //middle is the same as the name
  TString end = Form("%s%s",plotname.Data(),")");

  c0->Print(start.Data(),"pdf");
  c1->Print(plotname.Data(),"pdf"); 
  c2->Print(plotname.Data(),"pdf");
  c3->Print(plotname.Data(),"pdf");
  c4->Print(end.Data(),"pdf");

  //Write everything to output file
  fout->Write();
   
  // Send time efficiency report to console
  cout << "CPU time elapsed = " << watch->CpuTime() << " s = " << watch->CpuTime()/60.0 << " min. Real time = " << watch->RealTime() << " s = " << watch->RealTime()/60.0 << " min." << endl;

}//end Main
