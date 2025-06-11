#ifndef PHYSICS_CONSTANTS_H
#define PHYSICS_CONSTANTS_H

//Author: Ezekiel Wertz
//physics constants that we need for analysis
//Via pdg and pdg.lbl.gov/2020/reviews/rpp2020-rev-monte-carlo-numbering.pdf
namespace physics_constants{

  //math
  static const double PI = TMath::Pi();

  //photon
  static const double c = 299792458; // m/s
  static const int mc_idx_c = 22; //mc index

  //electron
  static const double q_e = 1.602176634E-19; // charge in Coulomb
  static const double M_e = 0.510998950E-03; // GeV/c^2
  static const int mc_idx_e = 11; //mc index

  //proton
  static const double M_p = 0.93827208816; // GeV/c^2
  static const int mc_idx_p = 2212; //mc index
  static const double mu_p = 2.79284735; //proton magnetic moment

  //neutron
  static const double M_n = 0.93956542052; // GeV/c^2
  static const int mc_idx_n = 2112; //mc index
  static const double mu_n = -1.9130427; ////Neutron anomalous magnetic moment

  //charged pion
  static const double M_pi = 0.13957; // GeV/c^2

  //atomic


  //other
  static const double alpha = 0.00729927; //unitless, fine structure constant
  
  
}
#endif
