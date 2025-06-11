#include "../include/physics.h"
#include "../include/physics_constants.h"

//Author: Ezekiel Wertz
//Companion implementation for the phyiscs functions

namespace physics{

 //Define HCal X-axis based on kinematic information
 TVector3 getHCal_xaxis(){
 //This is by convention
 TVector3 hcal_xaxis(0,-1,0);
 
 return hcal_xaxis;
 }

 //Define HCal Y-axis based on kinematic information
 TVector3 getHCal_yaxis(TVector3 hcal_xaxis, TVector3 hcal_zaxis){
 //Since we define hcal_xaxis and hcal_zaxis we require a perpendicular unit vector
 TVector3 hcal_yaxis = hcal_zaxis.Cross( hcal_xaxis ).Unit();

 return hcal_yaxis;
 }

 //Define HCal Z-axis based on kinematic information
 TVector3 getHCal_zaxis(double hcaltheta){
 //Define the zaxis so it is perpendicular to the face of HCal?
 TVector3 hcal_zaxis (sin(-hcaltheta),0,cos(-hcaltheta));
 
 return hcal_zaxis;
 }
 
 //Define the origin for HCal, part of defining the HCal coordinate system
 TVector3 getHCal_origin(double hcaldist, double hcal_offset, TVector3 hcal_xaxis, TVector3 hcal_zaxis){
 TVector3 hcal_origin = hcaldist *hcal_zaxis +hcal_offset*hcal_xaxis;

 return hcal_origin;
 }

 //Try to calculate the effects of the SBS magnetic field.
 double getBdL(int sbs_field,TString Kin,TString pass){
 double BdL = (exp_constants::getMaxSBSField(Kin,pass))*(exp_constants::sbsdipolegap)*(sbs_field/100.0);
 //cout << BdL << " " << exp_constants::getMaxSBSField(Kin) << " " << exp_constants::sbsdipolegap << " " << sbs_field << endl;
 return BdL;
 }

 //Get expected proton deflection, somewhat crude model
 //thetabend = 0.3*BdL/p
 double get_protonDeflection(double BdL, double p, double hcaldist, double sbsdist){
 double proton_deflection = tan(0.3*(BdL/p)) * (hcaldist - (sbsdist + (exp_constants::sbsdipolegap/2.0)));
 //cout << hcaldist << " " << sbsdist << " " << exp_constants::sbsdipolegap << endl;
 return proton_deflection;
 }

 //Calculate the outgoing energy based on losses in the target
 //Currently does not account for the aluminum target cell, requires modification
 double getEloss_outgoing(double bbtheta,TString target){
 double Eloss_outgoing;
        //LH2 
	if(target == "LH2"){
 	//Eloss_outgoing = (exp_constants::celldiameter/2.0)/sin(bbtheta)*(exp_constants::lh2_rho_tgt)*(exp_constants::lh2_dEdx); //Approximately 1 MeV, one could correct more using the raster position
	Eloss_outgoing = 0;
	//LD2
	}else if(target == "LD2"){
	//Eloss_outgoing = (exp_constants::celldiameter/2.0)/sin(bbtheta)*(exp_constants::ld2_rho_tgt)*(exp_constants::ld2_dEdx);
	Eloss_outgoing = 0;
	}else if(target == "Dummy"){
	//Not sure if this is correct kind of just using the others as an example and assuming its similar but with Aluminum for the windows.
	Eloss_outgoing = 0;
	//Eloss_outgoing = (exp_constants::celldiameter/2.0)/sin(bbtheta)*(exp_constants::rho_Al)*(exp_constants::dEdx_Al);
	}else{
	//give an error statement for and don't calculate the Eloss_outgoing.
	cout << "Warning: Target " << target << " is not handle by energy loss function! Defaulting to zero." << endl;
	Eloss_outgoing = 0;
	}

 return Eloss_outgoing;
 }

 //Calculate energy loss in target using vertex position
 double getEloss(double vz, TString target){
 double Eloss;

 	//LH2
 	if(target == "LH2"){
	//Eloss = ((exp_constants::l_tgt/2.0))*(exp_constants::lh2_rho_tgt)*(exp_constants::lh2_dEdx) + (exp_constants::lh2_uwallthick)*(exp_constants::rho_Al)*(exp_constants::dEdx_Al);
	Eloss = 0;
	//LD2
	}else if(target == "LD2"){
	//Eloss = ((exp_constants::l_tgt/2.0))*(exp_constants::ld2_rho_tgt)*(exp_constants::ld2_dEdx) + (exp_constants::ld2_uwallthick)*(exp_constants::rho_Al)*(exp_constants::dEdx_Al);
	Eloss = 0;
	//Not sure if this is correct kind of just using the others as an example and assuming its similar but with Aluminum for the windows.I
	}else if(target == "Dummy"){
	//Eloss = (exp_constants::ld2_uwallthick)*(exp_constants::rho_Al)*(exp_constants::dEdx_Al);
	Eloss = 0;
	}else{
        //give an error statement for and don't calculate the Eloss.
         cout << "Warning: Target " << target << " is not handle by energy loss function! Defaulting to zero." << endl;
	Eloss = 0;
	}
 return Eloss;	
 }

 //Correct beam energy for energy loss in target
 double getEcorr(double Ebeam, double Eloss){
 double Ecorr = Ebeam - Eloss; 
 return Ecorr;
 }
 
 //Make a vector based on vertex information. Right now this is only in vz. Could change later and include vx, and vy. Instead of assuming 0
 TVector3 getVertex(double vz){
 TVector3 vertex( 0.0, 0.0, vz );

 return vertex;
 }

 //reconstructed momentum, corrected for mean energy loss. Still need to include losses from Al shielding or target windows later
 double getp_recon_corr(double tr_p, double Eloss_outgoing){
 double pcorr = tr_p + Eloss_outgoing;
 return pcorr;
 }

 //four momentum vector for electron beam with correted Energy value
 TLorentzVector getpBeam(double Ecorr){
 TLorentzVector Pbeam(0.0,0.0,Ecorr,Ecorr);
 return Pbeam;
 }
  
 //four momentum for scattered electron based on reconstruction
 TLorentzVector getp_eprime(double tr_px, double tr_py, double tr_pz, double tr_p, double pcorr){
 TLorentzVector p_eprime(pcorr*(tr_px/tr_p),pcorr*(tr_py/tr_p),pcorr*(tr_pz/tr_p),pcorr);
 return p_eprime;
 }

 //four vector for target, currently implement for LD2 and LH2. Not sure if I should somehow implement just neutron
 TLorentzVector getp_targ(TString target){
 TLorentzVector ptarg;
 
 	//LH2 or protons
 	if(target == "LH2" || target == "p" || target == "H2"){
	ptarg.SetPxPyPzE(0.0,0.0,0.0,physics_constants::M_p);
	//LD2
	//For Quasi-elastic scattering dummy should be like deuterium
	}else if(target == "LD2" || target == "np"|| target == "Dummy" || target == "He3"){
	ptarg.SetPxPyPzE(0.0,0.0,0.0,0.5*(physics_constants::M_p+physics_constants::M_n));
	//just neutrons
	}else if(target == "n"){
	ptarg.SetPxPyPzE(0.0,0.0,0.0,physics_constants::M_n);
	}else{
	//give an error
	cout << "Error: Target " << target << " is not handled by this function! No value given." << endl;
	}
  return ptarg;
 }


 //four vector, virtual photon momentum or momentum transferred to the scattered nucleon
 TLorentzVector getq(TLorentzVector pbeam, TLorentzVector p_eprime){
 TLorentzVector q = pbeam - p_eprime;
 return q;
 }

 //Theta for scattered electron using reconstructed track momentum
 double get_etheta(TLorentzVector p_eprime){
 double etheta = acos(p_eprime.Pz()/p_eprime.E());
 return etheta;
 }

 //Phi for scattered electron using reconstructed track momentum
 double get_ephi(TLorentzVector p_eprime){
 double ephi = atan2(p_eprime.Py(),p_eprime.Px());
 /*cout << p_eprime.Py() << endl;
 cout << p_eprime.Px() << endl;
 cout << ephi << endl;*/
 return ephi;
 }

 //central momentum reconstructed from track angles and beam energy
 double get_pcentral(TLorentzVector pbeam,double etheta,TString target){
 double pcentral;
 double ebeam = pbeam.E();
	//LH2 or proton
	if(target == "LH2" || target == "p" || target == "H2"){
	pcentral = ebeam/(1.0 + (ebeam/physics_constants::M_p)*(1.0 - cos(etheta))); 
        //LD2 
        //For Quasi-elastic scattering dummy would be like deuterium
	}else if(target == "LD2" || target == "np"|| target == "Dummy" || target == "He3"){
	double Nmass = 0.5*(physics_constants::M_p+physics_constants::M_n);
        pcentral = ebeam/(1.0 + (ebeam/Nmass)*(1.0 - cos(etheta)));
        //neutron
	}else if(target == "n"){
	pcentral = ebeam/(1.0 + (ebeam/physics_constants::M_n)*(1.0 - cos(etheta)));
	}else{
        //give an error
	cout << "Error: Target " << target << " is not handled by this function! Defaulting to zero." << endl;
        pcentral = 0.0;
	}
 return pcentral;
 }

 //assume coplanarity, get the expected phi for the nucleon
 double get_phinucleon(double ephi,double PI){
 double phinucleon = ephi + PI;
 return phinucleon;
 }

 //Calculate Mott cross section for this event
 double getMott_CS(double alpha,double etheta,double pcorr, double Ecorr){
 double Mott_CS = (pow(alpha,2)*pow(cos(etheta/2.0),2)*pcorr)/(4.0*pow(Ecorr,3)*pow(sin(etheta/2.0),4));
 return Mott_CS;
 }

 //four momentum transferred squared, overload the function
 double getQ2(TLorentzVector q){
 double Q2 = (-q).M2(); 
 return Q2;
 }

 double getQ2(double ekineQ2){
 double Q2 = ekineQ2;
 return Q2;
 }

 double getQ2(TLorentzVector pbeam, TLorentzVector p_eprime, double etheta){
 double Q2 = 2.0*(pbeam.E())*(p_eprime.E())*(1.0 - cos(etheta));
 return Q2;
 }

 //four vector or unit vector, scattered nucleon momentum, overload the function
 TLorentzVector get_pN(TLorentzVector q,TLorentzVector p_targ){
 TLorentzVector p_N = q + p_targ;
 return p_N;
 }

 TLorentzVector get_pN(double p_N_exp,TVector3 p_Nhat, double nu, TLorentzVector p_targ){
 TLorentzVector p_N(p_N_exp*p_Nhat.X(),p_N_exp*p_Nhat.Y(),p_N_exp*p_Nhat.Z(),nu+p_targ.E());
 return p_N;
 }

 TVector3 get_pNhat(double theta_N_exp,double phi_N_exp){
 TVector3 p_Nhat(sin(theta_N_exp)*cos(phi_N_exp),sin(theta_N_exp)*sin(phi_N_exp),cos(theta_N_exp));
 return p_Nhat;
 }

 //energy transfer, overload the funtion
 double getnu(TLorentzVector q){
 double nu = q.E();
 return nu;
 }

 double getnu(double ekinenu){
 double nu = ekinenu;
 return nu;
 }

 double getnu(TLorentzVector pbeam, double pcentral){
 double nu = pbeam.E() - pcentral;
 return nu;
 }

 double getnu(TLorentzVector pbeam, TLorentzVector p_eprime){
 double nu = pbeam.E() - p_eprime.E();
 return nu;
 }

 //scattered nucleon expected momentum
 double get_pNexp(double nu,TString target){
 double p_N_exp;
 	 //LH2
	 if(target == "LH2"){
         p_N_exp = sqrt(pow(nu,2) + 2.0 * physics_constants::M_p * nu);
         //LD2
         //For Quasi-elastic scattering Dummy should be like deutreium
	 }else if(target == "LD2"||target == "Dummy"){
         double Nmass = 0.5*(physics_constants::M_p+physics_constants::M_n);
         p_N_exp = sqrt(pow(nu,2) + 2.0 * Nmass * nu);
         }else{
         //give an error
	 cout << "Error: Target " << target << " is not handled by this function! Defaulting to zero." << endl;
         p_N_exp = 0;
         }
 return p_N_exp;
 }

 //scattered nucleon expected angle
 double get_thetaNexp(TLorentzVector pbeam,double pcentral,double etheta,double p_N_exp){
 double theta_N_exp = acos((pbeam.E()-pcentral*cos(etheta))/p_N_exp);
 return theta_N_exp;
 }

 //Invariant Mass Squared, overloaded function
 double getW2(TLorentzVector p_N){
 double W2 = p_N.M2();
 return W2;
 }

 double getW2(double ekine_W2){
 double W2 = ekine_W2;
 return W2;
 }

 double getW2(TLorentzVector pbeam,TLorentzVector p_eprime, double Q2, TString target){
 double W2;
	//LH2
	if(target == "LH2" ){
	W2 = pow(physics_constants::M_p,2)+2.0*physics_constants::M_p*(pbeam.E() - p_eprime.E())-Q2;
	//LD2
	//For Quasi-elastic scattering dummy should be like deuterium
	}else if(target == "LD2"|| target == "Dummy"){
        double Nmass = 0.5*(physics_constants::M_p+physics_constants::M_n);
        W2 = pow(Nmass,2)+2.0*Nmass*(pbeam.E() - p_eprime.E())-Q2;
        }else{
        //give an error	
	cout << "Error: Target " << target << " is not handled by this function! Defaulting to zero." << endl;
        W2 = 0.0;
        }
 return W2;
 }

 //get the ray from Hall origin onto the face of hcal where the nucleon hit. This defines the intersection point of the nucleon with HCal
 TVector3 get_hcalintersect(TVector3 vertex,TVector3 hcal_origin,TVector3 hcal_zaxis,TVector3 p_Nhat ){
 // intersection of a ray with a plane in hcal coordinates
 double sintersect = ((hcal_origin - vertex).Dot(hcal_zaxis)) / (p_Nhat.Dot(hcal_zaxis));
 // ray from Hall origin onto the face of hcal where the nucleon hit
 TVector3 hcal_intersect = vertex + sintersect * p_Nhat;
 return hcal_intersect;
 }

 //gets expected location of scattered nucleon assuming straight line projections from BB track, x-direction
 double get_xhcalexpect(TVector3 hcal_intersect,TVector3 hcal_origin,TVector3 hcal_xaxis){
 double xhcal_expect = (hcal_intersect - hcal_origin).Dot(hcal_xaxis); 
 return xhcal_expect;
 } 

 //gets expected location of scattered nucleon assuming straight line projections from BB track, y-direction
 double get_yhcalexpect(TVector3 hcal_intersect,TVector3 hcal_origin,TVector3 hcal_yaxis){
 double yhcal_expect = (hcal_intersect - hcal_origin).Dot(hcal_yaxis);
 return yhcal_expect;
 }

 //intime cluster selection analysis,intime algorithm
 int cluster_intime_select(int num_hcal_clusid,double hcal_clus_atime[],double atime_sh,double hcal_clus_e[],double coin_mean,double coin_sig_fac,double coin_profile_sigma,double hcalemin){
 double maxE = 0.0;
 int intime_idx = -1;
 //loop through all clusters and select without HCal position information
 	for(int c = 0; c<num_hcal_clusid; c++){
	
	//calculate hcal physics quantities per cluster
	double atime = hcal_clus_atime[c];
        double atime_diff = atime - atime_sh; //Assuming best shower time on primary cluster
        double clus_energy = hcal_clus_e[c];
	//cout << c << " " << atime << " " << atime_sh << " " <<clus_energy << endl;
	//use hcal atime till after pass 2, wide cut around 5 sigma
	bool passCoin = abs(atime_diff - coin_mean) < coin_sig_fac*coin_profile_sigma;
	//cout << atime_diff << " " << coin_mean << " " << coin_sig_fac << " " << coin_profile_sigma << endl;
	//cout << c << " " << passCoin << " " << passE << endl;
	//in-time algorithm with new cluster, sort 
                if(passCoin){
		//if intime is passed see if the current cluster has the highest energy
		bool new_maxE = maxE < clus_energy;
			//if the current cluster has the higher energy update the information
			if(new_maxE){
			maxE = clus_energy;
			intime_idx = c;
			}
		}//end conditional
	}//end for loop over clusters
 return intime_idx;
 }

 //sort cluster to get highest energy cluster. This should just be a double check
 int cluster_HighEnergy(int num_hcal_clusid,double hcal_clus_e[]){
 double maxE = 0.0;
 int energy_idx = -1;
 	//loop through all the cluster
 	for(int c = 0; c<num_hcal_clusid; c++){
	double clus_energy = hcal_clus_e[c];

	bool new_maxE = maxE < clus_energy;

		//if the current cluster has the higher energy update the information
		if(new_maxE){
                maxE = clus_energy;
                energy_idx = c;
                }//end conditional	
	}//end loop

	//if(energy_idx > 0){
	//cout << "Index: " << energy_idx << " Value: " << maxE << endl;
	//}

 return energy_idx;
 }

 //sort clusters to get best intime indices from clone cluster, part 2 of intime algorithm
 int cluster_intime_findIdx(int num_hcal_clusid, vector<double> clone_cluster_intime){
 int intime_idx = -1;
 double intime = 0;

 	for(int d = 0;d < num_hcal_clusid; d++){

		if(clone_cluster_intime[d] > intime){
		intime = clone_cluster_intime[d];
		intime_idx = d;
		}//end conditional	
	}//end for loop

 return intime_idx;
 }

 //define dx
 double get_dx(double xhcal,double xhcal_expect){
 double dx = xhcal - xhcal_expect;
 return dx;
 }

 //define dy
 double get_dy(double yhcal,double yhcal_expect){
 double dy = yhcal - yhcal_expect;
 return dy;
 }

 //calculate the final MC weight best on generation iformation from simc
 double getMCWeight(double mc_weight,double luminosity,double genvol,int Ntried){

 double final_mc_weight = (mc_weight*luminosity*genvol)/Ntried;

 return final_mc_weight;

 }

 //calculate tau based on Q2 and particle type
 double get_tau(double Q2,TString target){
 double tau; 
 
 	//LH2 or protons
	if(target == "LH2" || target == "p" || target == "H2"){
	tau = Q2/(4.0*pow(physics_constants::M_p,2));
 	//LD2
 	//For Quasi-elastic scattering dummy should be like deuterium
 	}else if(target == "LD2" || target == "np"|| target == "Dummy" || target == "He3"){
	tau = Q2/(4.0*pow(0.5*(physics_constants::M_p + physics_constants::M_n),2));
 	//just neutrons
 	}else if(target == "n"){
	tau = Q2/(4.0*pow(physics_constants::M_n,2));
	}else{
        //give an error
	cout << "Error: Target " << target << " is not handled by this function! No value given." << endl;
        }
 return tau;
 }

 //calculate polarization of the virtual photon, epsilon
 double get_epsilon(double tau, double etheta){
 double epsilon = 1.0 / (1.0 + 2.0 * (1.0 + tau) * pow(tan(etheta/2.0),2));
 return epsilon;
 }

 //calculate the relative efficiency correction factor based on the given TH2D efficiency map and the position info from MC truth info
 double get_HCalEffCorr(TH2D* effMap, double x_mctrue, double y_mctrue, double accep_avg_eff,double pro_def, int fnucl){

 double x;

 //Get the bin content of the efficiency map at the given x and y from the MC truth info. That's the numerator
 //protons
 if(fnucl == 1){
 //account for the proton deflection
 x = x_mctrue - pro_def;
 //neutrons
 }else if(fnucl == 0){
 x = x_mctrue;
 }else{
 //we should never get here
 cout << "HCal Eff corr, we should never get here!" << endl;
 }

 
 //Convert the position info to bins
 
 int HCalxbin = effMap->GetYaxis()->FindBin(x);
 int HCalybin = effMap->GetXaxis()->FindBin(y_mctrue);

 //cout << "xbin: " << HCalxbin << " ybin: " << HCalybin << endl;

 int HCalxbin_cent = effMap->GetYaxis()->GetBinCenter(HCalxbin);
 int HCalybin_cent = effMap->GetXaxis()->GetBinCenter(HCalybin);

 
 double numerator = effMap->GetBinContent(HCalybin,HCalxbin);
 double denominator = accep_avg_eff;
 
 //cout << effMap->GetBinContent(HCalybin,HCalxbin) << endl;
 //cout << effMap->GetBinError(HCalybin,HCalxbin) << endl;

 double c = numerator/denominator;
 
 double bin_err = effMap->GetBinError(HCalybin,HCalxbin);

 //Low efficiency or empty region in the efficiency map
 //This decision is sort of empircial and limiting the bin_err essentially is a way to not apply a correction to bins with low statistics. Instead of directly cutting on the bin value which is the relative efficiency
 if(c == 0.0 || bin_err > 0.075){
 //if(c == 0.0 ){
 c =1.0;
 }

 //cout << c << endl;
 
 return c;
 }


}//end namepspace
