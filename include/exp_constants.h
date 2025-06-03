#ifndef EXP_CONSTANTS_H
#define EXP_CONSTANTS_H


//Author: Ezekiel Wertz
//constants specific to the GMn/nTPE experiment
namespace exp_constants{

  //Data
  static const int maxchan = 1000;
  static const int maxclus = 100;
  static const int minsamp = 0;
  static const int maxsamp = 40;     //total number of hcal ADC bins, GMn
  static const int maxsamps = 73000; //hcal chan * hcal max samp + 1000

  // Detectors
  // HCAL 
  // Dimensions
  static const int hcalchan = 288;
  static const int hcalcol = 12;
  static const int hcalrow = 24;
  static const double hcal_hrange = 1.85928;    //m, total range in horizontal direction of HCal (end-to-end)
  static const double hcal_vrange = 3.81;       //m, total range in vertical direction of HCal (end-to-end)
  static const double hcaladc_binw = 4.;        //ns, width of each ADC bin
  static const double hcalblk_w_p0 = 0.15;      //m, width of a HCAL block, use for pass0
  static const double hcalblk_h_p0 = 0.15;      //m, height of a HCAL block, use for pass0
  static const double hcalblk_w = 0.1524;       //m, width of a HCAL block, corrected
  static const double hcalblk_h = 0.1524;       //m, height of a HCAL block, corrected

  ////PASS2////
  static const double hcalblk_div_h = 0.15494;  //m, horizontal center-to-center dist.
  static const double hcalblk_div_v = 0.15875;  //m, vertical center-to-center dist.
  static const double hcalblk_div_hyp = 0.22183;//m, division corner-to-corner dist.
  static const double hcalblk_gap_h = 0.00254;  //m, horiz. gap bet. two blocks
  static const double hcalblk_gap_v = 0.00635;  //m, vert. gap bet. two blocks
 
  // Positions (mc)
  ////PASS2////
  static const double hcalposXi_mc = -2.27563;    //m, center of top row blocks (from SBS-replay/DB_MC)
  static const double hcalposXf_mc = 1.37562;     //m, center of bottom row blocks (from SBS-replay/DB_MC)
  static const double hcalposYi_mc = -0.85217;    //m, center of right most blocks (from SBS-replay/DB_MC)
  static const double hcalposYf_mc = 0.85217;     //m, center of left most blocks (from SBS-replay/DB_MC)
  
  // Pass0/1 (no block spacing)
  static const double hcalposXi_p0 = -2.16014;  //m, distance from beam center to top of HCal w/75cm offset
  static const double hcalposXf_p0 = 1.43826;   //m, distance from beam center to bottom of HCal w/75cm offset
  static const double hcalposYi_p0 = -0.9;      //m, distance from beam center to opposite-beam side of HCal
  static const double hcalposYf_p0 = 0.9;       //m, distance from beam center to beam side of HCal

  //BBCal
  static const int maxBBCalShChan = 189; // Total BBCal Shower Channels
  static const int maxBBCalShRows = 27;
  static const int maxBBCalShCols = 7;
  static const int maxBBCalPSChan = 52; // Total BBCal Preshower Channels
  static const int maxBBCalPSRows = 26;
  static const int maxBBCalPSCols = 2;

  //Beamline
  static const int chargeConvert = 3318; // See D.Flay Doc DB sbs.jlab.org/DocDB/0001/000164/002/dflay_bcm-ana-update_02-21-22.pdf p.8
  static const int clockActual = 103700; // Needed to convert the 104kHz clock to the actual counting rate
 
  //SBS Magnet
  static const double sbsdipolegap = 48.0*2.54/100.0; //about 1.22 m

  //////Static Target/Scattering Chamber Parameters
  // target
  static const double l_tgt = 60.0; // Length of the target (cm)
  static const double celldiameter = 1.82562; //cm, this is to properly cancel units, found using a picture of one of the target cells

  // H2
  static const double lh2_rho_tgt = 0.0723; // Density of target (g/cc)
  static const double lh2_cthick = 0.02;       //cm, target cell thickness
  static const double lh2_uwallthick = 0.0145; //cm, upstream wall thickness
  static const double lh2_dwallthick = 0.015;  //cm, downstream wall thickness
  static const double lh2_dEdx = 0.00574; //According to NIST ESTAR, the collisional stopping power of hydrogen is about 5.74 MeV*cm2/g at 2 GeV energy  
  static const double h2pressure = 1.013; //10 atm convert J/cm^3
  static const double h2temp = 303; //Room temperature

  // He3
  static const double he3pressure = 1.013;  //10 atm convert J/cm^3 - This should be different for each target...
  static const double he3temp = 305; // This seems wrong, He3 target was heated? Not 2deg above Room temp?

  //shielding
  static const double rho_Al = 2.7; // Density of aluminum windows (g/cc)
  static const double dEdx_Al = 0.0021; //According to NIST ESTAR, the collisional stopping power of Aluminum is about 2.1 MeV*cm2/g between 1-4 GeV
  static const double Alshieldthick = 2.54/8.0; //= 1/8 inch * 2.54 cm/inch 

  //Make functions to get hcal v offset. since it is pass dependent
  const double getHCalOffset( TString pass);
  
  const double getMaxSBSField(TString Kin,TString pass);

}
#endif
