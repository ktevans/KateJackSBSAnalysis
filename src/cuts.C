//cuts.C
//Author: Ezekiel Wertz
//The companion implementation for the header file
//A file to handle implementations for the cuts we need or care about for analysis

#include "../include/cuts.h"

namespace cuts {

//function to define HCal physical position boundaries
vector<double> hcal_Position_data(TString pass){
vector<double> hcalpos;

double hcal_Xi;
double hcal_Xf;
double hcal_Yi;
double hcal_Yf;

	if(pass == "pass0"){
	//slightly different calibration, need different parameters
	hcal_Xi = exp_constants::hcalposXi_p0;    
        hcal_Xf = exp_constants::hcalposXf_p0;
        hcal_Yi = exp_constants::hcalposYi_p0;
        hcal_Yf = exp_constants::hcalposYf_p0;
	//better calibration, differerent parameters
	}else if(pass == "pass1" || pass == "pass2"){
        hcal_Xi = exp_constants::hcalposXi_mc;      
        hcal_Xf = exp_constants::hcalposXf_mc;
        hcal_Yi = exp_constants::hcalposYi_mc;
        hcal_Yf = exp_constants::hcalposYf_mc;
        //catch all case in the event we get something not a pass we can handle
	}else{
        //put an error message here
	cout << "Error HCal Position  can't handle the pass: " << pass << "! Do something to fix it." << endl;
        }// end conditional
//store the values in a vector to use the information later
hcalpos.push_back(hcal_Xi);
hcalpos.push_back(hcal_Xf);
hcalpos.push_back(hcal_Yi);
hcalpos.push_back(hcal_Yf);

return hcalpos;
}

//function to define HCal physical position boundaries
vector<double> hcal_Position_MC(){
vector<double> hcalpos;

double hcal_Xi = exp_constants::hcalposXi_mc;
double hcal_Xf = exp_constants::hcalposXf_mc;
double hcal_Yi = exp_constants::hcalposYi_mc;
double hcal_Yf = exp_constants::hcalposYf_mc;

hcalpos.push_back(hcal_Xi);
hcalpos.push_back(hcal_Xf);
hcalpos.push_back(hcal_Yi);
hcalpos.push_back(hcal_Yf);

return hcalpos;
}


//functions to handle defining HCal active area
vector<double> hcal_ActiveArea_data(int num_blk_x, int num_blk_y, TString pass){

vector<double> hcalaa;
double hcalaa_Xi; 
double hcalaa_Xf;
double hcalaa_Yi;
double hcalaa_Yf;

	if(pass == "pass0"){
	//slightly different calibration, need different parameters
        hcalaa_Xi = exp_constants::hcalposXi_p0 + num_blk_x * exp_constants::hcalblk_h_p0;    
	hcalaa_Xf = exp_constants::hcalposXf_p0 - num_blk_x * exp_constants::hcalblk_h_p0;
	hcalaa_Yi = exp_constants::hcalposYi_p0 + num_blk_y * exp_constants::hcalblk_w_p0;
	hcalaa_Yf = exp_constants::hcalposYf_p0 - num_blk_y * exp_constants::hcalblk_w_p0;
	//better calibration, differerent parameters
	}else if(pass == "pass1" || pass == "pass2"){
	hcalaa_Xi = exp_constants::hcalposXi_mc + num_blk_x * exp_constants::hcalblk_div_v;      
	hcalaa_Xf = exp_constants::hcalposXf_mc - num_blk_x * exp_constants::hcalblk_div_v;
	hcalaa_Yi = exp_constants::hcalposYi_mc + num_blk_y * exp_constants::hcalblk_div_h;
	hcalaa_Yf = exp_constants::hcalposYf_mc - num_blk_y * exp_constants::hcalblk_div_h;
	//catch all case in the event we get something not a pass we can handle
	}else{
	//put an error message here
	cout << "Error HCal Active Area can't handle the pass: " << pass << "! Do something to fix it." << endl; 	
	}// end conditional

//store the values in a vector to use the information later
hcalaa.push_back(hcalaa_Xi);
hcalaa.push_back(hcalaa_Xf);
hcalaa.push_back(hcalaa_Yi);
hcalaa.push_back(hcalaa_Yf);

return hcalaa;
}//end function

vector<double> hcal_ActiveArea_MC(int num_blk_x, int num_blk_y){
vector<double> hcalaa;
double hcalaa_Xi = exp_constants::hcalposXi_mc + num_blk_x * exp_constants::hcalblk_div_v; 
double hcalaa_Xf = exp_constants::hcalposXf_mc - num_blk_x * exp_constants::hcalblk_div_v;
double hcalaa_Yi = exp_constants::hcalposYi_mc + num_blk_y * exp_constants::hcalblk_div_h;
double hcalaa_Yf = exp_constants::hcalposYf_mc - num_blk_y * exp_constants::hcalblk_div_h;

//store the values in a vector to use the information later
hcalaa.push_back(hcalaa_Xi);
hcalaa.push_back(hcalaa_Xf);
hcalaa.push_back(hcalaa_Yi);
hcalaa.push_back(hcalaa_Yf);

return hcalaa;
}//end function

//check position per event against hcal active area. Is True if detected on active area
bool hcalaa_ON (double xhcal, double yhcal, vector<double> hcalaa){

double hcalx_top = hcalaa[0];
double hcalx_bot = hcalaa[1];
double hcaly_left = hcalaa[2];
double hcaly_right = hcalaa[3];

//actually make the boolean value. Could be either true or false
bool detect = (yhcal > hcaly_left) && (yhcal < hcaly_right) && (xhcal > hcalx_top) && (xhcal < hcalx_bot); 

return detect;
}//end function

//function to define HCal fiducial area
vector<double> hcalfid(double dxsig_p, double dxsig_n, double dysig, vector<double> hcalaa,double num_sig_x, double num_sig_y ){
vector<double> hcalfid;

double hcalx_top = hcalaa[0] + num_sig_x * dxsig_p; //-x margin (relevant for proton)
double hcalx_bot = hcalaa[1] - num_sig_x * dxsig_n; //+x margin (relevant for neutron)
double hcaly_left = hcalaa[2] + num_sig_y * dysig; // -y margin
double hcaly_right = hcalaa[3] - num_sig_y * dysig; // +y margin

hcalfid.push_back(hcalx_top);
hcalfid.push_back(hcalx_bot);
hcalfid.push_back(hcaly_left);
hcalfid.push_back(hcaly_right);

return hcalfid;
}//end function

//function determine if event is inside the fiducial region. That is equally count neutrons and protons
bool hcalfid_IN(double xhcal_expect, double yhcal_expect, double dx_pn, vector<double> hcalfid){

double hcalx_top = hcalfid[0]; //-x fid
double hcalx_bot = hcalfid[1]; //+x fid
double hcaly_left = hcalfid[2]; //-y fid
double hcaly_right = hcalfid[3]; //+y fid

//for proton hypothesis
double xhcal_exp_pro = xhcal_expect - dx_pn; //define the proton expect postion based on the difference between the proton and neutron peaks

bool inFid = (yhcal_expect > hcaly_left) && (yhcal_expect < hcaly_right) &&  //check that dy is the same for protons and neutrons
	     (xhcal_expect > hcalx_top) && (xhcal_expect < hcalx_bot) &&  //check that expected position if it were a neutron
	     (xhcal_exp_pro > hcalx_top) && (xhcal_exp_pro < hcalx_bot);  //check that expected position if it were a proton

return inFid;
}//end function

//Funtion that defines good W2 elastic cut
bool goodW2(double W2, double W2_low, double W2_high){
bool goodW2 = (W2 >= W2_low) && (W2 <= W2_high);
return goodW2;
}

//Function used in defining globabl cut
bool failedGlobal(TTreeFormula *GlobalCut){
bool failedglobal = GlobalCut->EvalInstance(0) == 0;
return failedglobal;
}

//Function used to define good coin cut
bool passCoin(double coin_bestclus,double coin_mean, double coin_sig_fac, double coin_sigma){
bool passcoin = abs(coin_bestclus - coin_mean) <= coin_sig_fac*coin_sigma;
return passcoin;
}

//Function to define a good dy cut
bool good_dy(double dy_bestclus,double dyO_p, double dysig_cut_fac, double dysig_p){
bool good_dy = abs(dy_bestclus-dyO_p) <= dysig_cut_fac*dysig_p; 
return good_dy;
}

//Function to check if above min HCal E value
bool passHCalE(double hcal_e,double hcalemin){
bool passE = hcal_e >= hcalemin;
return passE;
}

//Function to check if number of cluster in HCal is above a min
bool passHCal_NClus(double nclus_hcal,int hcalnclusmin){
bool passNclus = nclus_hcal > hcalnclusmin;
return passNclus;
}


//A function that will determine how many sigma from the active area the expected value is. This function also incorporates information on if it is on the wrong side of the boundary, allowing negative values of sigma. Considers the Fiducial x-direction sigma factor.
double calculate_nsigma_fid_x(double xhcal_expect, double dxsig_p, double dxsig_n, double dx_pn,vector<double> hcalaa){
	//Active area top and bottom
	double hcalx_top = hcalaa[0];
	double hcalx_bot = hcalaa[1];

	double xhcal_exp_pro = xhcal_expect - dx_pn; //define the proton expect postion based on the difference between the proton and neutron peaks
	//Determine the distances to the top acceptance area boundary for both neutrons and protons
	double distance_to_top_x_aa_n = xhcal_expect - hcalx_top; // For neutron, this is because the top of HCal is in the -x direction
	double distance_to_top_x_aa_p = xhcal_exp_pro - hcalx_top; // For proton, this is because the top of HCal is in the -x direction

	//Now determine the number of sigma these values are from the acceptance boundary. Divide by dxsig_p because its the top
	double nsigma_top_n = distance_to_top_x_aa_n/dxsig_p;
	double nsigma_top_p = distance_to_top_x_aa_p/dxsig_p;

	//Determine the distances to the bottom acceptance area boundary for both neutrons and protons
	double distance_to_bot_x_aa_n = hcalx_bot - xhcal_expect; // For neutron, this is because the bot of HCal is in the +x direction
        double distance_to_bot_x_aa_p = hcalx_bot - xhcal_exp_pro; // For proton, this is because the bot of HCal is in the +x direction

	//Now determine the number of sigma these values are from the acceptance boundary. Divide by dxsig_n because its the bot
        double nsigma_bot_n = distance_to_bot_x_aa_n/dxsig_n;
        double nsigma_bot_p = distance_to_bot_x_aa_p/dxsig_n;

	//return the minimum of all the sigma values as that is the number of sigma we want
	return std::min({nsigma_bot_n,nsigma_bot_p,nsigma_top_n,nsigma_top_p});
}

//A function that will determine how many sigma from the active area the expected value is. This function also incorporates information on if it is on the wrong side of the boundary, allowing negative values of sigma. Considers the Fiducial y-direction sigma factor.
double calculate_nsigma_fid_y(double yhcal_expect, double dysig, vector<double> hcalaa){
	//Active Area left and right
	double hcaly_left = hcalaa[2];
	double hcaly_right = hcalaa[3];

	//Determine the distances to the left (negative) acceptance area boundary
	double distance_to_left_y_aa = yhcal_expect - hcaly_left; //The left of HCal should be from the -y direction
	//Now determine the number of sigma these values are from the acceptance boundary.
	double nsigma_left = distance_to_left_y_aa / dysig;

	//Determine the distances to the right (positive) acceptance area boundary
	double distance_to_right_y_aa = hcaly_right - yhcal_expect; //The right of HCal should be from the +y direction
	//Now determine the number of sigma these values are from the acceptance boundary.
	double nsigma_right = distance_to_right_y_aa / dysig;

	//return the minimum of all the sigma values as that is the number of sigma we want
	return std::min(nsigma_left,nsigma_right);
}

//nsigma fid check
bool passNsigFid(double nsigx_fid,double nsigy_fid){
bool passNsigFid = nsigx_fid > 1 && nsigy_fid >1;

return passNsigFid;
}

//Spot cut on HCal dx and dy to select protons or neutrons
bool passHCal_Spot(double dx, double dy, double dxO, double dyO, double dxsigma, double dysigma, double radius){
bool HCal_spot = pow((dx-dxO)/dxsigma,2.0) + pow((dy-dyO)/dysigma,2.0) <= pow(radius,2.0);
return HCal_spot;
}

}//end namespace
