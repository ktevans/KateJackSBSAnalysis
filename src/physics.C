#ifndef PHYSICS_H
#define PHYSICS_H

#include "TLorentzVector.h"
#include "TVector3.h"
//Author Ezekiel Wertz
//A namespace to hold functions that calculate or update relevant physics analysis variables or quantities


namespace physics{

 //Define HCal X-axis based on kinematic information
 TVector3 getHCal_xaxis();

 //Define HCal Y-axis based on kinematic information
 TVector3 getHCal_yaxis(TVector3 hcal_xaxis, TVector3 hcal_yaxis);

 //Define HCal Z-axis based on kinematic information
 TVector3 getHCal_zaxis(double hcaltheta);

 //Define the origin for HCal, part of defining the HCal coordinate system
 TVector3 getHCal_origin(double hcaldist, double hcal_offset, TVector3 hcal_xaxis, TVector3 hcal_zxaxis);

 //Try to calculate the effects of the SBS magnetic field
 double getBdL(int sbs_field,TString Kin,TString pass);

 //Get expected proton deflection, somewhat crude model
 //thetabend = 0.3*BdL/p
 double get_protonDeflection(double BdL, double p, double hcaldist, double sbsdist);

 //Calculate the outgoing energy based on losses in the target
 //Currently does not account for the aluminum target cell, requires modification
 double getEloss_outgoing(double bbtheta,TString target);

 //Calculate the energy loss in target using vertex position
 double getEloss(double vz, TString target);

 //Correct beam energy for energy loss in target
 double getEcorr(double Ebeam, double Eloss);
 
 //Make a vector based on vertex information. Right now this is only in vz. Could change later and include vx, and vy. Instead of assuming 0
 TVector3 getVertex(double vz); 

 //reconstructed momentum, corrected for mean energy loss. Still need to include losses from Al shielding or target windows later
 double getp_recon_corr(double tr_p, double Eloss_outgoing); 
 
 //four momentum vector for electron beam with correted Energy value
 TLorentzVector getpBeam(double Ecorr);

 //four momentum for scattered electron based on reconstruction
 TLorentzVector getp_eprime(double tr_px, double tr_py, double tr_pz, double tr_p, double pcorr);

 //four vector for target, currently implement for LD2 and LH2. Implement for just protons and neutrons as well
 TLorentzVector getp_targ(TString target);

 //four vector, virtual photon momentum or momentum transferred to the scattered nucleon
 TLorentzVector getq(TLorentzVector pbeam, TLorentzVector p_eprime);

 //Theta for scattered electron using reconstructed track momentum
 double get_etheta(TLorentzVector p_eprime);

 //Phi for scattered electron using reconstructed track momentum
 double get_ephi(TLorentzVector p_eprime);

 //central momentum reconstructed from track angles and beam energy
 double get_pcentral(TLorentzVector pbeam,double etheta,TString target);

 //assume coplanarity, get the expected phi for the nucleon
 double get_phinucleon(double ephi,double PI);

 //Calculate Mott cross section for this event
 double getMott_CS(double alpha, double etheta, double pcorr, double Ecorr);

 //four momentum transferred squared, overload the function
 double getQ2(TLorentzVector q);

 double getQ2(double ekineQ2);

 double getQ2(TLorentzVector pbeam, TLorentzVector p_eprime, double etheta);

 //four vector or unit vector, scattered nucleon momentum, overload the function
 TLorentzVector get_pN(TLorentzVector q,TLorentzVector p_targ);
 
 TLorentzVector get_pN(double p_N_exp,TVector3 p_Nhat, double nu, TLorentzVector p_targ);

 TVector3 get_pNhat(double theta_N_exp,double phi_N_exp);

 //energy transfer, overload the funtion
 double getnu(TLorentzVector q);

 double getnu(double ekinenu);

 double getnu(TLorentzVector pbeam, double pcentral);

 double getnu(TLorentzVector pbeam, TLorentzVector p_eprime);

 //scattered nucleon expected momentum
 double get_pNexp(double nu,TString target);

 //scattered nucleon expected angle
 double get_thetaNexp(TLorentzVector pbeam,double pcentral,double etheta,double p_N_exp);

 //Invariant Mass Squared, overloaded function
 double getW2(TLorentzVector p_N);

 double getW2(double ekine_W2); 
 
 double getW2(TLorentzVector pbeam,TLorentzVector p_eprime, double Q2, TString target);

 //get the ray from Hall origin onto the face of hcal where the nucleon hit. This defines the intersection point of the nucleon with HCal
 TVector3 get_hcalintersect(TVector3 vertex,TVector3 hcal_origin,TVector3 hcal_zaxis,TVector3 p_Nhat );

 //gets expected location of scattered nucleon assuming straight line projections from BB track, x-direction
 double get_xhcalexpect(TVector3 hcal_intersect,TVector3 hcal_origin,TVector3 hcal_xaxis);

 //gets expected location of scattered nucleon assuming straight line projections from BB track, y-direction
 double get_yhcalexpect(TVector3 hcal_intersect,TVector3 hcal_origin,TVector3 hcal_yaxis);

 //intime cluster selection analysis, part 1 of intime algorithm
 int cluster_intime_select(int num_hcal_clusid,double hcal_clus_atime[],double atime_sh,double hcal_clus_e[],double coin_mean,double coin_sig_fac,double coin_profile_sigma,double hcalemin);

 //sort clusters to get best intime indices from clone cluster, part 2 of intime algorithm
 int cluster_intime_findIdx(int num_hcal_clusid, vector<double> clone_cluster_intime);

 //define dx
 double get_dx(double xhcal,double xhcal_expect);

 //define dy
 double get_dy(double yhcal,double yhcal_expect);

 //sort cluster to get highest energy cluster. This should just be a double check
 int cluster_HighEnergy(int num_hcal_clusid,double hcal_clus_e[]);

 //calculate the final MC weight best on generation iformation from simc
 double getMCWeight(double mc_weight,double luminosity,double genvol,int Ntried);

 //calculate tau based on Q2 and particle type
 double get_tau(double Q2,TString target);

 //calculate polarization of the virtual photon, epsilon
 double get_epsilon(double tau, double etheta);

 //calculate the relative efficiency correction factor based on the given TH2D efficiency map and the position info from MC truth info
 double get_HCalEffCorr(TH2D* effMap, double x_mctrue, double y_mctrue, double accep_avg_eff,double pro_def, int fnucl);

}//end namespace

#endif
