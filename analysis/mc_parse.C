//Author Ezekiel Wertz
//03/15/2024
//Purpose: Parsing simulation information from a simc generator by nucleon, build histrograms for further analysis

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
#include "../src/utility.C"
#include "../src/exp_constants.C"
#include "../src/kinematic_obj.C"
#include "../src/data_object.C"
#include "../src/cuts.C"
#include "../src/physics.C"
#include "../src/parse_config.C"

//Main
void mc_parse(const char *setup_file_name){

//Define a clock to check macro processing time
TStopwatch *watch = new TStopwatch();
watch->Start( kTRUE );

//parse object to get in the information that The One Config file has and is manipulated
parse_config mainConfig(setup_file_name);
//mainConfig.printMCYields();

//store all the parameters from the mainConfig file into local variables. So we don't have to keep recalling them
TString rootfile_dir = mainConfig.getRootFileDir();
TString histfile_dir = mainConfig.getHistFileDir();
TCut globalcut = mainConfig.getGlobalCut();
TString exp = mainConfig.getExp();
TString kin = mainConfig.getKin();
TString kinematic_file = mainConfig.getKinFileName();
int sbs_field = mainConfig.getSBSField();
TString partial_name_p = mainConfig.getPartialNameP();
TString partial_name_n = mainConfig.getPartialNameN();
TString pass = mainConfig.getPass();
TString target = mainConfig.getTarg();
bool sync_jobs = mainConfig.get_syncJobs();
bool mc_override = mainConfig.get_MCOverride();
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
double dxsig_n_fac = mainConfig.get_dxSignFac();
double dxsig_p_fac = mainConfig.get_dxSigpFac();
double dysig_n_fac = mainConfig.get_dySignFac();
double dysig_p_fac = mainConfig.get_dySigpFac();
double W2fitmax = mainConfig.getW2FitMax();
double binfac = mainConfig.getBinFac();
double hbinfac = mainConfig.getHBinFac();
int e_method = mainConfig.get_emethod();
int maxtracks = mainConfig.getMAXNTRACKS();
double hcalemin = mainConfig.getHCaleMin();
double dysig_cut = mainConfig.get_dySigCut();
double dysig_cut_fac = mainConfig.get_dySigCutFac();
int hcalnclusmin = mainConfig.get_HCalNclusMin();
double hcal_offset = exp_constants::getHCalOffset(pass);

//config info for efficiency map
//bool HCal_Eff_map_flag = mainConfig.get_HCalEffMap();
//double HCal_accep_avg_eff = mainConfig.get_HCalAccepAvgEff();
//TString HCal_Eff_map_file = mainConfig.get_HCalEffMapFile();

//double hcalfit_low = exp_constants::hcalposXi_mc; //lower fit/bin limit for hcal dx plots. 
//double hcalfit_high = exp_constants::hcalposXf_mc; //higher fit/bin limit for hcal dx plots.
//double hcal_fitrange = exp_constants::hcal_vrange; //Full range of hcal dx plots

//store all important kinematic info in local variables
kinematic_obj myKin(kinematic_file,kin);
double Ebeam = myKin.getBeamEnergy();
double hcaldist = myKin.getHCalDist();
double sbsdist = myKin.getSBSDist();
double bbtheta = myKin.getBBAngle_Rad();
double hcaltheta = myKin.getHCalAngle_Rad();
double p_nuc_centr = myKin.getNucleonP();

//setup hcal physical bounds that match database
vector<double> hcalpos = cuts::hcal_Position_MC();

//setup hcal active area with bounds that match database
vector<double> hcalaa = cuts::hcal_ActiveArea_MC(1,1);

//setup input file for the efficiency map
//TFile *eff_map_file = new TFile(HCal_Eff_map_file.Data());

//setup output file
TString outfile = utility::makeOutputFileName_MCParse(exp,kin,sbs_field,target); 
TFile *fout = new TFile(outfile,"RECREATE");

//need to find the MC root and hist files 
//proton
vector<string> rootFileNames_p;
vector<string> histFileNames_p;

//neutron
vector<string> rootFileNames_n;
vector<string> histFileNames_n;

//metadata proton
vector<pair<string,vector<float>>> metadata_p;

//metadata neutron
vector<pair<string,vector<float>>> metadata_n;

//Need to allow the MC parser to handle LD2, He3, and LH2 for MC file handling

if(target == "LD2" || target == "He3")
{
	//find information from CSV file and root file paths from generic sim replay output supported by jlab-HPC
	//for proton
	utility::SyncFilesCSV(histfile_dir,rootfile_dir,partial_name_p,rootFileNames_p,metadata_p);
	//for neutron
	utility::SyncFilesCSV(histfile_dir,rootfile_dir,partial_name_n,rootFileNames_n,metadata_n);    
}
else if(target == "LH2" || target == "H2")
{
	//find information from CSV file and root file paths from generic sim replay output supported by jlab-HPC
        //for proton
        utility::SyncFilesCSV(histfile_dir,rootfile_dir,partial_name_p,rootFileNames_p,metadata_p);
}
else 
{
	cout << "Error: During MC file import. The target type: " << target << " is not supported. Resolve issue!" << endl;
}
//double check that we have the same number of files again
int pFiles = 0;
int nFiles = 0;

//Need to allow the MC parser to handle both LD2 and LH2 for MC file handling

if(target == "LD2" || target == "He3")
{
	//Continuing to make sure there exist both protons and neutrons. Make sure the job has a both a proton and neutron file. If not get rid of it
	if(sync_jobs)
	{
		//Function that searchs, finds, and then removes any files that are not in both proton and neutron vectors
		//cout << metadata_p.size() << " " << metadata_n.size() << endl;
		utility::syncJobNumbers(metadata_p,metadata_n);
		pFiles = metadata_p.size();
		nFiles = metadata_n.size();
	}	
   
	cout << endl << "Number of proton files " << pFiles << " , Number of neutron files " << nFiles << endl;
	if(pFiles != nFiles)
	{
		cout << endl << "Warning: Number proton and neutron root files loaded does not match. Investigate!" << endl << endl;
		if(sync_jobs)
		{
			return;
		}
	}
}
else if(target == "LH2" || target == "H2")
{
	//Don't need to sync jobs. But do need to make sure the different replay types populate the number of files
	pFiles = metadata_p.size();
}
else 
{
	cout << "Error: During MC file import. The target type: " << target << " is not supported. Resolve issue!" << endl;
}

double dx_pn = mainConfig.get_dxpn();

  //create output tree
  TTree *Parse = new TTree("Parse","Analysis Tree");

  //uncut output tree variables

  //calculated and CODA variables
  double dx_out, dy_out, xexp_out, yexp_out, W2_out, Q2_out, nu_out, tau_out, epsilon_out, pcorr_out, mott_out;
  Parse->Branch("dx", &dx_out, "dx/D");
  Parse->Branch("dy", &dy_out, "dy/D");
  Parse->Branch("sbs.hcal.x_exp", &xexp_out, "sbs.hcal.x_exp/D");
  Parse->Branch("sbs.hcal.y_exp", &yexp_out, "sbs.hcal.y_exp/D");
  Parse->Branch("e.kine.W2", &W2_out, "e.kine.W2/D");
  Parse->Branch("e.kine.Q2", &Q2_out, "e.kine.Q2/D")
  Parse->Branch("e.kine.nu", &nu_out, "e.kine.nu/D");
  Parse->Branch("e.kine.tau", &tau_out, "e.kine.tau/D");
  Parse->Branch("e.kine.epsilon", &epsilon_out, "e.kine.epsilon/D");
  Parse->Branch("pcorr", &pcorr_out, "pcorr/D");
  Parse->Branch("e.kine.mott", &mott_out, "e.kine.mott/D");

  //hcal variables
  double sbs_hcal_x_out, sbs_hcal_y_out, sbs_hcal_e_out, sbs_hcal_atimeblk_out, sbs_hcal_clus_blk_row_out, sbs_hcal_clus_blk_col_out, sbs_hcal_nblk_out;
  Parse->Branch("sbs.hcal.x", &sbs_hcal_x_out, "sbs.hcal.x/D");
  Parse->Branch("sbs.hcal.y", &sbs_hcal_y_out, "sbs.hcal.y/D");
  Parse->Branch("sbs.hcal.e", &sbs_hcal_e_out, "sbs.hcal.e/D");
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
  
  double coin_mean_out;
  double coin_sigma_out;

  //kinematic variables
  int file_out, mag_out;
  Parse->Branch("filenum", &file_out, "filenum/I");
  Parse->Branch("sbs.field", &mag_out, "sbs.field/I");

  double proton_deflection_out, p_N_out, p_central_out, MC_p_n_out, MC_p_e_out, MC_vz_out;
  Parse->Branch("proton_deflection", &proton_deflection_out, "proton_deflection/D");
  Parse->Branch("p_N", &p_N_out, "p_N/D");
  Parse->Branch("p_central", &p_central_out, "p_centrial/D");
  Parse->Branch("MC_p_n", &MC_p_n_out, "MC_p_n/D");
  Parse->Branch("MC_p_e", &MC_p_e_out, "MC_p_e/D");
  Parse->Branch("MC_vz", &MC_vz_out, "MC_vz/D");

  //cut variables
  int passGlobal_out, HCalON_out;
  Parse->Branch("passGlobal", &passGlobal_out,"passGlobal/I");
  Parse->Branch("HCalON", &HCalON_out,"HCalON/I");

  double hcal_sh_atime_diff_out;
  Parse->Branch("adc.coin", &hcal_sh_atime_diff_out, "adc.coin/D");
  
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
  //double uncorr_mc_weight_out;
  //double corr_mc_weight_out;
  //double x_mctrue_n_out;
  //double y_mctrue_n_out;
  //double p_mctrue_n_out;
  //double q_mctrue_n_out;
  //double theta_mctrue_n_out;
  //double phi_mctrue_n_out;
  //double pcentral_mctrue_n_out;
  //double phi_N_exp_out;
  //double Q2_mctrue_n_out;
  //double nu_mctrue_n_out;
  //double theta_N_exp_out;
  //double p_Nhat_mctrue_n_out;
  //double p_N_mctrue_n_out;
  //double W2_mctrue_n_out;
  //double p_Nhat_out;
  //double tau_mctrue_n_out;
  //double epsilon_mctrue_n_out;
  //double hcal_intersect_out;
  //double hcal_intersect_mctrue_n_out;
  //double dx_mctrue_n_out;
  //double dy_mctrue_n_out;
  //int passW2_out;
  //int passCoin_out;
  //int passFid_out;
  //int is_proton_out;
  //int is_neutron_out;
  //double corr_fac_out;

  //logistical information
  TString nuc = "none"; 
  int num_nuc = 0; 
  int num_files = 0;
  int is_proton = 0; //false
  int is_neutron = 0; //false

  //Conditional to handle the num_nuc depending on target
  if(target == "LD2" || target == "He3")
  {
	num_nuc = 2;	
  }
  else if(target == "LH2" || target == "H2")
  {
	num_nuc = 1;
  }
  else
  {
	//stays default zero
	cout << "Error: During event loop setup. The target type: " << target << " is not supported. Resolve issue!" << endl;	
  }

  TChain *C = nullptr;

  //Main loop over nucleons (r==0 proton, r==1 neutron)
  for(int r=0; r<num_nuc; r++)
  {

  	//update current nucleon
  	if(r==0)
 	{
		nuc = "p";
		is_proton = 1;
		is_neutron = 0;
		num_files = pFiles;
	}
	else if(r==1)
	{
		nuc = "n";
		is_proton = 0;
		is_neutron = 1;
		num_files = nFiles;
	}

	//loop over simulation files for given iteration
	for(int f=0; f<num_files; ++f)
	{
		string root_file;
		string hist_file;
		int Ntried,job_number,Nthrown;
                double luminosity,genvol,ebeam,charge,seed,weight_gt_max,obs_max_weight;
		bool using_rej_samp = 0;
		double max_weight = 0;
		
		//get information for MC weights
		if(r==0)
		{
			if(mc_override)
			{
				Ntried = mainConfig.getNTriedOverride();
                		luminosity = mainConfig.getLumiOverride();
                		genvol = mainConfig.getVolOverride();
			}
			else if(replay_type == "jlab-HPC")
			{
				root_file = metadata_p[f].first;
				//Verify with CSV file structure
				job_number = metadata_p[f].second[0];
				Nthrown = metadata_p[f].second[1];
				Ntried = metadata_p[f].second[2];		
				genvol = metadata_p[f].second[3];
				luminosity = metadata_p[f].second[4];
				ebeam = metadata_p[f].second[5];
				charge = metadata_p[f].second[6];
				seed = metadata_p[f].second[7];
				//Make sure we have conditional for old vs new files
				if(metadata_p[f].second.size() > 8)
				{
					using_rej_samp = metadata_p[f].second[8];
					max_weight = metadata_p[f].second[9];
					weight_gt_max = metadata_p[f].second[10];
					obs_max_weight = metadata_p[f].second[11];
				}
			}
			else
			{
				cerr << "No MC weight information provided. Invesitgate this problem! No values assigned." << endl;
			}
		}
		else if(r==1)
		{
			if(mc_override)
			{
                        	Ntried = mainConfig.getNTriedOverride();
                        	luminosity = mainConfig.getLumiOverride();
                        	genvol = mainConfig.getVolOverride();
                        }
			else if(replay_type == "jlab-HPC")
			{
                        	root_file = metadata_n[f].first;
                        	//Verify with CSV file structure
				job_number = metadata_n[f].second[0];
                        	Nthrown = metadata_n[f].second[1];
                        	Ntried = metadata_n[f].second[2];
                        	genvol = metadata_n[f].second[3];
                        	luminosity = metadata_n[f].second[4];
                        	ebeam = metadata_n[f].second[5];
                        	charge = metadata_n[f].second[6];
                        	seed = metadata_n[f].second[7];
				//Make sure we have conditional for old vs new files
				if(metadata_n[f].second.size() > 8)
				{
                        		using_rej_samp = metadata_n[f].second[8];
                        		max_weight = metadata_n[f].second[9];
                       			weight_gt_max = metadata_n[f].second[10];
                                	obs_max_weight = metadata_n[f].second[11];
                                }
                        }
			else
			{
                        	cerr << "No MC weight information provided. Invesitgate this problem! No values assigned." << endl;
                        }
		}
		cout << "For this " << nuc << " file: Ntried = " << Ntried << " luminosity = " << luminosity << " genvol = " << genvol << " max weight " << max_weight << endl;

		//add the file to the TChain
		C = new TChain("T");
		C->Add(root_file.c_str());

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


		//Monteo Carlo Variables
		double mc_weight;
		C->SetBranchStatus("MC.simc_Weight",1);

		C->SetBranchAddress("MC.simc_Weight", &mc_weight);

		//BBCal shower
		double atime_sh, e_sh, nclus_sh, nblk_sh, BB_E_over_p;

  		C->SetBranchStatus( "bb.sh.atimeblk", 1 );
  		C->SetBranchStatus( "bb.sh.e", 1 );
  		C->SetBranchStatus( "bb.sh.nclus", 1 );
  		C->SetBranchStatus( "bb.sh.nblk", 1 );
  		C->SetBranchStatus( "bb.etot_over_p", 1);

  		C->SetBranchAddress("bb.sh.atimeblk", &atime_sh);
  		C->SetBranchAddress("bb.sh.e", &e_sh);
  		C->SetBranchAddress("bb.sh.nclus", &nclus_sh);
  		C->SetBranchAddress("bb.sh.nblk", &nblk_sh);
  		C->SetBranchAddress( "bb.etot_over_p", &BB_E_over_p );

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

		//Branches for MC truth info related to Efficiency map
		double MC_true_p_n, MC_true_px_n, MC_true_py_n, MC_true_pz_n,MC_true_p_e, MC_true_px_e, MC_true_py_e, MC_true_pz_e,MC_true_vz, MC_true_Ebeam;
	        double mc_fnucl; //1 = proton, 0 = neutron
		C->SetBranchStatus("MC.simc_p_n",1);
		C->SetBranchStatus("MC.simc_px_n",1);
		C->SetBranchStatus("MC.simc_py_n",1);
		C->SetBranchStatus("MC.simc_pz_n",1);
		C->SetBranchStatus("MC.mc_fnucl",1);
		C->SetBranchStatus("MC.simc_p_e",1);
                C->SetBranchStatus("MC.simc_px_e",1);
                C->SetBranchStatus("MC.simc_py_e",1);
                C->SetBranchStatus("MC.simc_pz_e",1);
		C->SetBranchStatus("MC.simc_vz",1);
		C->SetBranchStatus("MC.simc_Ebeam",1);

		C->SetBranchAddress("MC.simc_p_n",&MC_true_p_n);
		C->SetBranchAddress("MC.simc_px_n",&MC_true_px_n);
		C->SetBranchAddress("MC.simc_py_n",&MC_true_py_n);
		C->SetBranchAddress("MC.simc_pz_n",&MC_true_pz_n);
		C->SetBranchAddress("MC.mc_fnucl",&mc_fnucl);
		C->SetBranchAddress("MC.simc_p_e",&MC_true_p_e);
                C->SetBranchAddress("MC.simc_px_e",&MC_true_px_e);
                C->SetBranchAddress("MC.simc_py_e",&MC_true_py_e);
                C->SetBranchAddress("MC.simc_pz_e",&MC_true_pz_e);
		C->SetBranchAddress("MC.simc_vz",&MC_true_vz);
		C->SetBranchAddress("MC.simc_Ebeam",&MC_true_Ebeam);


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
		while(C->GetEntry(nevent++))
		{
			//progress tracker
			cout << "Processing " << nuc  << " file " << f << "/" << num_files << " event " << nevent << "/" << nentries << "\r";
        		cout.flush();
			//single loop global cut
			currenttreenum = C->GetTreeNumber();
        		if( nevent == 1 || currenttreenum != treenum )
			{
        			treenum = currenttreenum;
        			GlobalCut->UpdateFormulaLeaves();
        		}

        		//Is true if failed global cut
			bool failglobal = cuts::failedGlobal(GlobalCut);
			///////////
			//Electron-arm physics calculations

			//BB E/p calculated. Should be the same as the tree variable. But it does not cause seg faults
        		double BB_E_over_p = (e_sh+e_ps)/tr_p[0];

			//use the ebeam from the MC file if possible, if not will default to value from kinematic info
			if((ebeam/1000) != -1 )
			{
				Ebeam = ebeam/1000;
			}
						
			//correct beam energy from vertex information
			double Eloss = physics::getEloss(tr_vz[0],target);
        		double Ecorr = physics::getEcorr(Ebeam,Eloss);
			
			//make the vertex, assuming only beam direction
			TVector3 vertex = physics::getVertex(tr_vz[0]);
			//reconstructed momentum from track momentum information, corrected for mean energy loss. Still need to include losses from Al shielding or target windows later
			//Vertex but for the MC true info
			TVector3 vertex_mctrue_n = physics::getVertex(MC_true_vz);
			
			//double pcorr = physics::getp_recon_corr(tr_p[0],Eloss_outgoing);
			double pcorr = tr_p[0];

			//four momentum vector for electron beam with correted Energy value
			TLorentzVector pbeam = physics::getpBeam(Ebeam);

			//four momentum for scattered electron based on reconstruction
			TLorentzVector p_eprime = physics::getp_eprime(tr_px[0],tr_py[0],tr_pz[0],tr_p[0],pcorr);

			//four momentum for efficiency correction from MC truth info
			TLorentzVector p_mctrue_n = physics::getp_eprime(MC_true_px_n,MC_true_py_n,MC_true_pz_n,MC_true_p_n,MC_true_p_n);

			//four vector for target
			TLorentzVector p_targ = physics::getp_targ(nuc);

			//four vector, virtual photon momentum or momentum transferred to the scattered nucleon
			TLorentzVector q = physics::getq(pbeam,p_eprime);
        		TVector3 q_vec = q.Vect();

			//q vector associated with the four momentum transfer for the MC truth info
			TLorentzVector q_mctrue_n = physics::getq(pbeam, p_mctrue_n);
			TVector3 q_mctrue_n_vec = q_mctrue_n.Vect(); 

			//Theta for scattered electron using reconstructed track momentum
			double etheta = physics::get_etheta(p_eprime);

			//Phi for scattered electron using reconstructed track momentum
			double ephi = physics::get_ephi(p_eprime);

			//Determine the theta and phi associated with the MC truth info of the nucleon. Use these directly instead of the projected values
			double theta_mctrue_n = physics::get_etheta(p_mctrue_n);
			double phi_mctrue_n = physics::get_ephi(p_mctrue_n);
			//Because of coordinate system
			phi_mctrue_n = physics_constants::PI/2.0 - phi_mctrue_n;

			//central momentum reconstructed from track angles and beam energy
			double pcentral = physics::get_pcentral(pbeam,etheta,nuc);

			//central momentum using the MC true info for the nucleon
			double pcentral_mctrue_n = physics::get_pcentral(pbeam,theta_mctrue_n,nuc);

			//assume coplanarity, get the expected phi for the nucleon
			double phi_N_exp = physics::get_phinucleon(ephi,physics_constants::PI);

			//Calculate Mott cross section for this event
			double Mott_CS = physics::getMott_CS(physics_constants::alpha,etheta,pcorr,Ebeam);

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

			//The same directly above variables, just implemented using the MC truth info for nucleon instead, part of efficiency correction. Here we are going to assume we are using the angles method which is method 3.
			//enable the same calculation but for the MC truth info, part of efficiency correction calculation
			double Q2_mctrue_n = physics::getQ2(pbeam,p_mctrue_n,theta_mctrue_n);
                        double nu_mctrue_n = physics::getnu(pbeam,pcentral_mctrue_n);
                        double p_N_exp_mctrue_n = physics::get_pNexp(nu_mctrue_n,target);
                        TVector3 p_Nhat_mctrue_n = physics::get_pNhat(theta_mctrue_n,phi_mctrue_n);
                        TLorentzVector p_N_mctrue_n = physics::get_pN(p_N_exp_mctrue_n,p_Nhat_mctrue_n,nu_mctrue_n,p_targ);
                        double W2_mctrue_n = physics::getW2(pbeam,p_mctrue_n,Q2_mctrue_n,target);
                        double tau_mctrue_n = physics::get_tau(Q2_mctrue_n,target);
                        double epsilon_mctrue_n = physics::get_epsilon(tau_mctrue_n,theta_mctrue_n);
	
			//conditional to determine remaining e-arm related calculations.
			if(e_method == 1)
			{
				//v1
				Q2 = physics::getQ2(q);
                		p_N = physics::get_pN(q,p_targ);
                		p_Nhat = p_N.Vect().Unit();
                		nu = physics::getnu(q);
                		W2 = physics::getW2(p_N);
                		tau = physics::get_tau(Q2,target);
                		epsilon = physics::get_epsilon(tau,etheta);			
			}
			else if(e_method == 2)
			{
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
			}
			else if(e_method == 3)
			{
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
			}
			else if(e_method == 4)
			{
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
			}
			else
			{
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
			TVector3 hcal_intersect = physics::get_hcalintersect(vertex,hcal_origin,hcal_zaxis,p_Nhat);

			//HCal intersect but from quantities derived from MC truth information
			TVector3 hcal_intersect_mctrue_n = physics::get_hcalintersect(vertex_mctrue_n,hcal_origin,hcal_zaxis,p_Nhat_mctrue_n);
				
			//gets expected location of scattered nucleon assuming straight line projections from BB track, x-direction
			double xhcal_expect = physics::get_xhcalexpect(hcal_intersect,hcal_origin,hcal_xaxis);

			//gets expected location of scattered nucleon assuming straight line projections from BB track, y-direction
			double yhcal_expect = physics::get_yhcalexpect(hcal_intersect,hcal_origin,hcal_yaxis);

			//Now calculate x_expect and y_expect with the MC derived truth info
			double xhcal_mctrue_n = physics::get_xhcalexpect(hcal_intersect_mctrue_n,hcal_origin,hcal_xaxis);
			//Because of coordinate system difference
			xhcal_mctrue_n = -xhcal_mctrue_n;
			double yhcal_mctrue_n = physics::get_yhcalexpect(hcal_intersect_mctrue_n,hcal_origin,hcal_yaxis);
				

			//Calculate expected proton deflection with a somewhat crude module
        		double BdL = physics::getBdL(sbs_field,kin,pass);
        		double proton_deflection = physics::get_protonDeflection(BdL,p_N.Vect().Mag(),hcaldist,sbsdist);
			double proton_deflection_mctrue_n = physics::get_protonDeflection(BdL,p_N_mctrue_n.Vect().Mag(),hcaldist,sbsdist);

			//supposedly the MC does timing poorly. So implementing an intime clustering algorithm for MC is not the best idea. Just do a highest energy cluster search to be sure.
			int energy_idx = physics::cluster_HighEnergy(num_hcal_clusid,hcal_clus_e);

			//Assume the highest energy cluster sort is sufficient to find the best cluster
			int clus_idx_best = energy_idx;
				

			//calculate important information
			double xhcal_bestclus = hcal_clus_x[clus_idx_best];
        		double yhcal_bestclus = hcal_clus_y[clus_idx_best];
        		double dx_bestclus = physics::get_dx(xhcal_bestclus,xhcal_expect);
        		double dy_bestclus = physics::get_dy(yhcal_bestclus,yhcal_expect);
        		double hcal_atime_bestclus = hcal_clus_atime[clus_idx_best];
        		double coin_bestclus = hcal_atime_bestclus - atime_sh;
        		double coin_pclus = hcal_clus_atime[0] - atime_sh;
        		double hcal_e_bestclus = hcal_clus_e[clus_idx_best];
			int hcal_nblk_bestclus = (int) hcal_clus_nblk[clus_idx_best];

			//Check on the MC truth info by calculating dx and dy
			double dx_mctrue_n = physics::get_dx(xhcal_bestclus,xhcal_mctrue_n);
                        double dy_mctrue_n = physics::get_dy(yhcal_bestclus,yhcal_mctrue_n);

			//calculate the number of sigma away from fiducial boundaries. Store info for later
        		double nsigx_fid = cuts::calculate_nsigma_fid_x(xhcal_expect,dxsig_p,dxsig_n,dx_pn,hcalaa);
        		double nsigy_fid = cuts::calculate_nsigma_fid_y(yhcal_expect,dysig_p,hcalaa);
				
			//setup booleans for cuts later. Save boolean values to tree

			//HCal active area
			bool hcalaa_ON = cuts::hcalaa_ON(xhcal_bestclus,yhcal_bestclus,hcalaa);

			//W2 elastic boolean
			bool goodW2 = cuts::goodW2(W2,W2_low,W2_high);

			//good dy boolean
			bool good_dy = cuts::good_dy(dy_bestclus,dyO_p,dysig_cut_fac,dysig_cut);

			//good coincidence time cut. MC timing is poor do not implement this cut.
			//bool passCoin = cuts::passCoin(coin_pclus,coin_mean,coin_sig_fac,coin_profile_sig);

			//good fiducial cut
			//bool passFid = cuts::hcalfid_IN(xhcal_expect,yhcal_expect,dx_pn,hcalfid);

			//pass HCal E
			bool passHCalE = cuts::passHCalE(hcal_e_bestclus,hcalemin);
				
			//pass HCal num clus
			bool passHCal_Nclus = cuts::passHCal_NClus(nclus_hcal,hcalnclusmin);

			//pass NSig Fid check
        		//bool passNSigFid = cuts::passNsigFid(nsigx_fid,nsigy_fid);

			//calculate the final corrected MC weight
			//handles both cases of using rejection sampling and not 
			double uncorr_mc_weight;
			double corr_mc_weight;
			if(using_rej_samp)
			{
				uncorr_mc_weight = physics::getMCWeight(max_weight,luminosity,genvol,Ntried);
			}
			else
			{
				uncorr_mc_weight = physics::getMCWeight(mc_weight,luminosity,genvol,Ntried);
			}
				
			//If a efficiency map correction is involved apply it to the final_mc_weight
				
			double corr_fac = 0.0;

			//Temp change to check what happens if we only altered just the proton files
			if(HCal_Eff_map_flag)
			{
				//if(HCal_Eff_map_flag && r==0){
				corr_fac = physics::get_HCalEffCorr(xy_expect_eff_map_clone,xhcal_mctrue_n,yhcal_mctrue_n,HCal_accep_avg_eff,proton_deflection_mctrue_n,mc_fnucl);
				corr_mc_weight = uncorr_mc_weight * corr_fac;
			}
			else
			{
				corr_mc_weight = uncorr_mc_weight;
			}

			//Fill analysis tree variables before making cuts
			dx_out = dx_bestclus;
        		dy_out = dy_bestclus;
        		xexp_out = xhcal_expect;
        		yexp_out = yhcal_expect;
        		xhcal_out = xhcal_bestclus;
        		yhcal_out = yhcal_bestclus;
        		W2_out = W2;
        		Q2_out = Q2;
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
        		//passFid_out = (int) passFid;
        		file_out = f;
        		mag_out = sbs_field;
			uncorr_mc_weight_out = uncorr_mc_weight;
			corr_mc_weight_out = corr_mc_weight;
			is_proton_out = is_proton;
			is_neutron_out = is_neutron;
			hcal_sh_atime_diff_out = coin_bestclus;
        		proton_deflection_out = proton_deflection;
      			p_N_out = p_N.Vect().Mag();
        		p_central_out = pcentral;
        		nblkHCAL_out = nblkHCAL;
        		rowblkHCAL_out = rowblkHCAL[clus_idx_best];
        		colblkHCAL_out = colblkHCAL[clus_idx_best];
			MC_p_n_out = MC_true_p_n;
			MC_p_e_out = MC_true_p_e;
			MC_vz_out = MC_true_vz;

			x_mctrue_n_out = xhcal_mctrue_n;
			y_mctrue_n_out = yhcal_mctrue_n;
			p_mctrue_n_out = p_mctrue_n.Vect().Mag();
			q_mctrue_n_out = q_mctrue_n.Vect().Mag();
			theta_mctrue_n_out = theta_mctrue_n;
			phi_N_exp_out = phi_N_exp;
			phi_mctrue_n_out = phi_mctrue_n;
			pcentral_mctrue_n_out = pcentral_mctrue_n;
			Q2_mctrue_n_out = Q2_mctrue_n;
			nu_mctrue_n_out = nu_mctrue_n;
			theta_N_exp_out = theta_N_exp;
			p_Nhat_mctrue_n_out = p_Nhat_mctrue_n.Mag2();
			p_N_mctrue_n_out = p_N_mctrue_n.Vect().Mag();
			W2_mctrue_n_out = W2_mctrue_n;
			p_Nhat_out = p_Nhat.Mag2();
			tau_mctrue_n_out = tau_mctrue_n;
			epsilon_mctrue_n_out = epsilon_mctrue_n;
			hcal_intersect_out = hcal_intersect.Mag2();
			hcal_intersect_mctrue_n_out = hcal_intersect_mctrue_n.Mag2();
			corr_fac_out = corr_fac;

			dx_mctrue_n_out = dx_mctrue_n;
			dy_mctrue_n_out = dy_mctrue_n;
				
			//Fill the analysis tree
			Parse->Fill();
		}//end event loop
		C->Reset();
	}//end loop simulation
  }//end nucleon for loop

  //Write everything to output file
  fout->Write();

// Send time efficiency report to console
cout << "CPU time elapsed = " << watch->CpuTime() << " s = " << watch->CpuTime()/60.0 << " min. Real time = " << watch->RealTime() << " s = " << watch->RealTime()/60.0 << " min." << endl;

}//end main
