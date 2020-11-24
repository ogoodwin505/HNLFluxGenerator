// #ifndef HNLFlux_funcs_cxx
#define HNLFlux_funcs_cxx


#include "HNLFlux_funcs.h"
#include <iostream>
#include <iomanip>
#include <string>

#include "TChain.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TRotation.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TGraph.h"

#include <cmath>
#include <vector>
#include <map>
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"




  static const int kpdg_nue_G       =   53;  // extended Geant 53
  static const int kpdg_nuebar_G   =  52;  // extended Geant 52
  static const int kpdg_numu_G      =   56;  // extended Geant 56
  static const int kpdg_numubar_G   =  55;  // extended Geant 55

  // static const int kpdg_muplus_G     =   5;  // Geant  5
  // static const int kpdg_muminus_G    =    6;  // Geant  6
  // static const int kpdg_pionplus_G   =   8;  // Geant  8
  // static const int kpdg_pionminus_G  =  9;  // Geant  9
  // static const int kpdg_k0long_G     =   10;  // Geant 10  ( K0=311, K0S=310 )
  // static const int kpdg_k0short_G    =   16;  // Geant 16
  // static const int kpdg_k0mix_G      =   311;
  // static const int kpdg_kaonplus_G   =   11;  // Geant 11
  // static const int kpdg_kaonminus_G  =  12;  // Geant 12
  // static const int kpdg_omegaminus_G =  24;  // Geant 24
  // static const int kpdg_omegaplus_G  = 32;  // Geant 32



  
  const int kpdg_nue       =   12;  // extended Geant 53
  const int kpdg_nuebar    =  -12;  // extended Geant 52
  const int kpdg_numu       =   14;  // extended Geant 53
  const int kpdg_numubar    =  -14;  // extended Geant 52

  const int kpdg_eplus     =    -11;  // Geant  5
  const int kpdg_eminus    =     11;  // Geant  6
  const int kpdg_muplus     =   -13;  // Geant  5
  const int kpdg_muminus    =    13;  // Geant  6
  const int kpdg_pionplus   =   211;  // Geant  8
  const int kpdg_pionminus  =  -211;  // Geant  9
  const int kpdg_pi0        =   111;  // Geant  9
  const int kpdg_k0long     =   130;  // Geant 10  ( K0=311, K0S=310 )
  const int kpdg_k0short    =   310;  // Geant 16 //not used none of the dk2nu files
  const int kpdg_k0mix      =   311;  //not used
  const int kpdg_kaonplus   =   321;  // Geant 11
  const int kpdg_kaonminus  =  -321;  // Geant 12
  const int kpdg_omegaminus =  3334;  // Geant 24 dont get used
  const int kpdg_omegaplus =  -3334;  // Geant 24 dont get used

  const int kpdg_muplus_G     =   5;  // Geant  5
  const int kpdg_muminus_G    =    6;  // Geant  6
  const int kpdg_pionplus_G   =   8;  // Geant  8
  const int kpdg_pionminus_G  =  9;  // Geant  9
  const int kpdg_k0long_G     =   10;  // Geant 10  ( K0=311, K0S=310 )
  const int kpdg_k0short_G    =   16;  // Geant 16
  const int kpdg_k0mix_G      =   311;
  const int kpdg_kaonplus_G   =   11;  // Geant 11
  const int kpdg_kaonminus_G  =  12;  // Geant 12
  const int kpdg_omegaminus_G =  24;  // Geant 24
  const int kpdg_omegaplus_G  = 32;  // Geant 32



  const double kPIMASS = 0.13957;
  const double kKMASS  = 0.49368;
  const double kK0MASS = 0.49767;
  const double kMUMASS = 0.105658389;
  const double kOMEGAMASS = 1.67245;


  bool debug1 = false;

  std::map<int, std::map<double, TGenPhaseSpace*>> decayers;


//___________________________________________________________________________
int Random_e_mu(Long64_t counter) {
  //This function is used to return e or mu like events half the time each
  TRandom3 *r = new TRandom3(counter);

  int output = 0;
  double e_nu = r->Uniform(0.0, 10.0);

  if (e_nu < 5.0) {
    output = 2;
  }
  else {
    output = 3;
  }

  delete r;

  return output;
}

//___________________________________________________________________________
TVector3 RandomInTPC(Long64_t counter, Parameter& par) {

  TRandom3 *r = new TRandom3(counter);




  TVector3 UpperBound = par.DetUpperBound();
  TVector3 LowerBound = par.DetLowerBound();

  double x = r->Uniform(LowerBound.X(), UpperBound.X());
  double y = r->Uniform(LowerBound.Y(), UpperBound.Y());
  double z = r->Uniform(LowerBound.Z(), UpperBound.Z());



  TVector3 det;
  det.SetXYZ(x, y, z);
  delete r;

  return det;
}



//___________________________________________________________________________
TVector3 FromDetToBeam( const TVector3& det, Parameter& par) {
  TVector3 beam;

  // TRotation R;
  // TVector3 newX(0.921228671,   0.00136256111, -0.389019125);
  // TVector3 newY(0.0226872648,  0.998103714,    0.0572211871);
  // TVector3 newZ(0.388359401,  -0.061539578,    0.919450845);
  //
  // R.RotateAxes(newX,newY,newZ);
  // if (debug1) {
  //     cout << "R_{beam to det} = " << endl;
  //     cout << " [ " << R.XX() << " " << R.XY() << " " << R.XZ() << " ] " << endl;
  //     cout << " [ " << R.YX() << " " << R.YY() << " " << R.YZ() << " ] " << endl;
  //     cout << " [ " << R.ZX() << " " << R.ZY() << " " << R.ZZ() << " ] " << endl;
  //     cout << endl;
  // }
  // R.Invert(); // R is now the inverse
  // if (debug1) {
  //     cout << "R_{det to beam} = " << endl;
  //     cout << " [ " << R.XX() << " " << R.XY() << " " << R.XZ() << " ] " << endl;
  //     cout << " [ " << R.YX() << " " << R.YY() << " " << R.YZ() << " ] " << endl;
  //     cout << " [ " << R.ZX() << " " << R.ZY() << " " << R.ZZ() << " ] " << endl;
  //     cout << endl;
  //     }
  //Now R allows to go from detector to beam coordinates.
  //NuMIDet is vector from NuMI target to uB detector (in beam coordinates)
  //
  //Shouldnt need rotation for bnb

  // TVector3 BNBDet (-124.325, 0.93, 46336.35); //cm
  beam = det + par.DetPos_vec();

  return beam;
}

//___________________________________________________________________________
double estimate_pots(int highest_potnum) {

  // Stolen: https://cdcvs.fnal.gov/redmine/projects/dk2nu/repository/show/trunk/dk2nu
  // looks like low counts are due to "evtno" not including
  // protons that miss the actual target (hit baffle, etc)
  // this will vary with beam conditions parameters
  // so we should round way up, those generating flugg files
  // aren't going to quantize less than 1000
  // though 500 would probably be okay, 100 would not.
  // also can't use _last_ potnum because muons decay->> don't
  // have theirs set
  //
  // Marco: Trying with 10000

  const Int_t    nquant = 10000; //1000; // 500;  // 100
  const Double_t rquant = nquant;

  Int_t estimate = (TMath::FloorNint((highest_potnum - 1) / rquant) + 1) * nquant;
  return estimate;
}

//___________________________________________________________________________
int calcEnuWgt(bsim::Dk2Nu* dk2nu, const TVector3& xyz,double& enu, double& wgt_xy) {

  const double kPIMASS = 0.13957;
  const double kKMASS  = 0.49368;
  const double kK0MASS = 0.49767;
  const double kMUMASS = 0.105658389;
  const double kOMEGAMASS = 1.67245;

  const double kRDET = 100.0;   // set to flux per 100 cm radius

  double xpos = xyz.X();
  double ypos = xyz.Y();
  double zpos = xyz.Z();

  enu    = 0.0;  // don't know what the final value is
  wgt_xy = 0.0;  // but set these in case we return early due to error


  double parent_mass = kPIMASS;
  switch ( dk2nu->decay.ptype ) {
  case kpdg_pionplus:
  case kpdg_pionminus:
    parent_mass = kPIMASS;
    break;
  case kpdg_kaonplus:
  case kpdg_kaonminus:
    parent_mass = kKMASS;
    break;
  case kpdg_k0long:
  case kpdg_k0short:
  case kpdg_k0mix:
    parent_mass = kK0MASS;
    break;
  case kpdg_muplus:
  case kpdg_muminus:
    parent_mass = kMUMASS;
    break;
  case kpdg_omegaminus:
  case kpdg_omegaplus:
    parent_mass = kOMEGAMASS;
    break;
  default:
    std::cerr << "bsim::calcEnuWgt unknown particle type " << dk2nu->decay.ptype
              << std::endl << std::flush;
    //assert(0);
    return 1;
  }
  double parentp2 = ( dk2nu->decay.pdpx * dk2nu->decay.pdpx +
                      dk2nu->decay.pdpy * dk2nu->decay.pdpy +
                      dk2nu->decay.pdpz * dk2nu->decay.pdpz );
  double parent_energy = TMath::Sqrt( parentp2 +
                                      parent_mass * parent_mass);
  double parentp = TMath::Sqrt( parentp2 );

  double gamma     = parent_energy / parent_mass;
  double gamma_sqr = gamma * gamma;
  double beta_mag  = TMath::Sqrt( ( gamma_sqr - 1.0 ) / gamma_sqr );

  // Get the neutrino energy in the parent decay CM
  double enuzr = dk2nu->decay.necm;
  // Get angle from parent line of flight to chosen point in beam frame
  double rad = TMath::Sqrt( (xpos - dk2nu->decay.vx) * (xpos - dk2nu->decay.vx) +
                            (ypos - dk2nu->decay.vy) * (ypos - dk2nu->decay.vy) +
                            (zpos - dk2nu->decay.vz) * (zpos - dk2nu->decay.vz) );

  double emrat = 1.0;
  double  costh_pardet = -999.;
  //double theta_pardet = -999.;
  // boost correction, but only if parent hasn't stopped
  if ( parentp > 0. ) {
    costh_pardet = ( dk2nu->decay.pdpx * (xpos - dk2nu->decay.vx) +
                     dk2nu->decay.pdpy * (ypos - dk2nu->decay.vy) +
                     dk2nu->decay.pdpz * (zpos - dk2nu->decay.vz) )
                   / ( parentp * rad);
    if ( costh_pardet >  1.0 ) costh_pardet =  1.0;
    if ( costh_pardet < -1.0 ) costh_pardet = -1.0;
    //theta_pardet = TMath::ACos(costh_pardet);
    // Weighted neutrino energy in beam, approx, good for small theta

    emrat = 1.0 / ( gamma * ( 1.0 - beta_mag * costh_pardet ));
  }
  enu = emrat * enuzr;  // the energy ... normally
  //std::cout << enu << endl;
  // Get solid angle/4pi for detector element
  //     // small angle approximation, fixed by Alex Radovic
  //         //SAA//  double sangdet = ( kRDET*kRDET /
  //             //SAA//                   ( (zpos-decay->Vz)*(zpos-decay->Vz) ) ) / 4.0;
  //

  double sanddetcomp = TMath::Sqrt( ( (xpos - dk2nu->decay.vx) * (xpos - dk2nu->decay.vx) ) +
                                    ( (ypos - dk2nu->decay.vy) * (ypos - dk2nu->decay.vy) ) +
                                    ( (zpos - dk2nu->decay.vz) * (zpos - dk2nu->decay.vz) )   );
  double sangdet = (1.0 - TMath::Cos(TMath::ATan( kRDET / sanddetcomp ))) / 2.0;

  // Weight for solid angle and lorentz boost
  wgt_xy = sangdet * ( emrat * emrat );  // ! the weight ... normally
  // Done for all except polarized muon decay
  //     // in which case need to modify weight
  //         // (must be done in double precision)

  if ( dk2nu->decay.ptype == kpdg_muplus || dk2nu->decay.ptype == kpdg_muminus) {
    double beta[3], p_dcm_nu[4], p_nu[3], p_pcm_mp[3], partial;

    // Boost neu neutrino to mu decay CM
    beta[0] = dk2nu->decay.pdpx / parent_energy;
    beta[1] = dk2nu->decay.pdpy / parent_energy;
    beta[2] = dk2nu->decay.pdpz / parent_energy;
    p_nu[0] = (xpos - dk2nu->decay.vx) * enu / rad;
    p_nu[1] = (ypos - dk2nu->decay.vy) * enu / rad;
    p_nu[2] = (zpos - dk2nu->decay.vz) * enu / rad;
    partial = gamma *
              (beta[0] * p_nu[0] + beta[1] * p_nu[1] + beta[2] * p_nu[2] );
    partial = enu - partial / (gamma + 1.0);
    // the following calculation is numerically imprecise
    // especially p_dcm_nu[2] leads to taking the difference of numbers
    // of order ~10's and getting results of order ~0.02's
    // for g3numi we're starting with floats (ie. good to ~1 part in 10^7)

    p_dcm_nu[0] = p_nu[0] - beta[0] * gamma * partial;
    p_dcm_nu[1] = p_nu[1] - beta[1] * gamma * partial;
    p_dcm_nu[2] = p_nu[2] - beta[2] * gamma * partial;
    p_dcm_nu[3] = TMath::Sqrt( p_dcm_nu[0] * p_dcm_nu[0] +
                               p_dcm_nu[1] * p_dcm_nu[1] +
                               p_dcm_nu[2] * p_dcm_nu[2] );

    // Boost parent of mu to mu production CM

    double particle_energy = dk2nu->decay.ppenergy;
    gamma = particle_energy / parent_mass;
    beta[0] = dk2nu->decay.ppdxdz * dk2nu->decay.pppz / particle_energy;
    beta[1] = dk2nu->decay.ppdydz * dk2nu->decay.pppz / particle_energy;
    beta[2] =                    dk2nu->decay.pppz / particle_energy;
    partial = gamma * ( beta[0] * dk2nu->decay.muparpx +
                        beta[1] * dk2nu->decay.muparpy +
                        beta[2] * dk2nu->decay.muparpz );
    partial = dk2nu->decay.mupare - partial / (gamma + 1.0);
    p_pcm_mp[0] = dk2nu->decay.muparpx - beta[0] * gamma * partial;
    p_pcm_mp[1] = dk2nu->decay.muparpy - beta[1] * gamma * partial;
    p_pcm_mp[2] = dk2nu->decay.muparpz - beta[2] * gamma * partial;
    double p_pcm = TMath::Sqrt ( p_pcm_mp[0] * p_pcm_mp[0] +
                                 p_pcm_mp[1] * p_pcm_mp[1] +
                                 p_pcm_mp[2] * p_pcm_mp[2] );

    const double eps = 1.0e-30;  // ? what value to use
    if ( p_pcm < eps || p_dcm_nu[3] < eps ) {
      return 3; // mu missing parent info?
    }
    // Calc new decay angle w.r.t. (anti)spin direction
    double costh = ( p_dcm_nu[0] * p_pcm_mp[0] +
                     p_dcm_nu[1] * p_pcm_mp[1] +
                     p_dcm_nu[2] * p_pcm_mp[2] ) /
                   ( p_dcm_nu[3] * p_pcm );
    if ( costh >  1.0 ) costh =  1.0;
    if ( costh < -1.0 ) costh = -1.0;
    // Calc relative weight due to angle difference
    double wgt_ratio = 0.0;
    switch ( dk2nu->decay.ntype ) {
    case kpdg_nue:
    case kpdg_nuebar:
      wgt_ratio = 1.0 - costh;
      break;
    case kpdg_numu:
    case kpdg_numubar:
    {
      double xnu = 2.0 * enuzr / kMUMASS;
      wgt_ratio = ( (3.0 - 2.0 * xnu )  - (1.0 - 2.0 * xnu) * costh ) / (3.0 - 2.0 * xnu);
      break;
    }
    default:
      return 2; // bad neutrino type
    }

    wgt_xy = wgt_xy * wgt_ratio;

  } // ptype is muon
  return 0;
}

//___________________________________________________________________________
double findPathLen(bsim::Dk2Nu* dk2nu, const TVector3& xyz) {
  const TVector3 pos(xyz.X() - dk2nu->decay.vx, xyz.Y() - dk2nu->decay.vy, xyz.Z() - dk2nu->decay.vz);
  double length = pos.Mag();
  return length;

}

//___________________________________________________________________________
double DalitzDensity(const std::map<int, double> pLambda_map, const std::map<int, double> pXi0_map, const TLorentzVector& p4_k, const TLorentzVector& p4_pi, const TLorentzVector& p4_lep, const TLorentzVector& p4_Hnu, double decay_type) {

  double E_pi = p4_pi.E();
  double E_lep = p4_lep.E();
  double E_Hnu = p4_Hnu.E();

  double massK = p4_k.M();
  double massPi = p4_pi.M();
  double massLep = p4_lep.M();
  double massHnu = p4_Hnu.M();

  double E_pi_max = 0.5 * (massK * massK + massPi * massPi - massLep * massLep - massHnu * massHnu) / massK;
  double E = E_pi_max - E_pi;
  double Q2 = massK * massK + massPi * massPi - 2. * massK * E_pi;

  const double pLambda = pLambda_map.at(decay_type);
  const double pXi0 = pXi0_map.at(decay_type);


  double F = 1. + pLambda * Q2 / massPi / massPi;
  //double Fmax = 1.;
  //  //if(pLambda > 0.)
  double Fmax = (1. + pLambda * (massK * massK / massPi / massPi + 1.));

  double Xi = pXi0 * (1. + pLambda * Q2 / massPi / massPi);

  double coeffA = massK * (2. * E_lep * E_Hnu - massK * E) + massLep * massLep * (E / 4. - E_Hnu);
  double coeffB = massLep * massLep * (E_Hnu - E / 2.);
  double coeffC = massLep * massLep * E / 4.;

  double RhoMax = (Fmax * Fmax) * (massK * massK * massK / 8.);
  double Rho = (F * F) * (coeffA + coeffB * Xi + coeffC * Xi * Xi);

  return Rho / RhoMax;
}


//___________________________________________________________________________
std::map<int, double> get_Q_vals(const std::map<int, std::vector<int>>& daughters, const std::map<int, int> mother_pdgs) {
  std::map<int, double> Q_vals;
  for (auto decay_def : daughters) {
    double Q = TDatabasePDG::Instance()->GetParticle(mother_pdgs.at(decay_def.first))->Mass();
    for (auto dpdg : decay_def.second) {
      Q -=  TDatabasePDG::Instance()->GetParticle(dpdg)->Mass();
    }
    Q_vals[decay_def.first] = Q;
  }
  return Q_vals;
}


//___________________________________________________________________________
TGenPhaseSpace* get_decayer(int decay_type, double heavy_mass, const std::map<int, std::vector<int>> daughter_pdgs, const std::map<int, int> mother_pdgs) {

  if (daughter_pdgs.find(decay_type) == daughter_pdgs.end()) {
    throw "Decayer failed";
  }
  std::map<double, TGenPhaseSpace*>& decay = decayers[decay_type];
  if (decay.find(heavy_mass) != decay.end()) {
    return decay[heavy_mass];
  }

  std::vector<double> daughter_masses = { heavy_mass };
  for (auto dpdg : daughter_pdgs.at(decay_type)) {
    daughter_masses.push_back(TDatabasePDG::Instance()->GetParticle(dpdg)->Mass());
  }

  TGenPhaseSpace *decayer = new TGenPhaseSpace;
  TLorentzVector dummy(0., 0., 0., TDatabasePDG::Instance()->GetParticle(mother_pdgs.at(decay_type))->Mass());
  decayer->SetDecay(dummy, daughter_masses.size(), daughter_masses.data());
  decay[heavy_mass] = decayer;
  return decayer;
}


//___________________________________________________________________________
TLorentzVector make_decay(const std::map<int, double> pLambda_map, const std::map<int, double> pXi0_map, TGenPhaseSpace* decayer, const TLorentzVector& parentP4,  int decay_type, double *Decay_weight) {
  const double wmax = decayer->GetWtMax();
  double weight = -1.;
  const unsigned int nmaxthrows = 100000;
  unsigned int nthrows = 0;

  while (gRandom->Uniform(wmax) > weight) {
    weight = decayer->Generate();

    // Kaon dalitz decay
    double ddw_weight = weight;
    if (decay_type <= 4 || decay_type == 6 || decay_type == 7 || decay_type == 9 || decay_type == 10 ) {
      double ddw = DalitzDensity(pLambda_map, pXi0_map, parentP4, *(decayer->GetDecay(1)), *(decayer->GetDecay(2)), *(decayer->GetDecay(0)), decay_type);

      if (gRandom->Uniform() > ddw) {
        weight = -1.; // will cause rethrow at while(gRandom->Uniform(wmax) > weight)
        ddw_weight *= ddw;
      }
    }

    // muon decay
    if (decay_type == 11 || decay_type == 12) {
      // TODO: Muon decay beta-style energy spectrum
    }

    // break infintie loop
    if (nthrows++ > nmaxthrows) {

      if (weight < 0.) {
        // kaon dalitz decay
        weight = ddw_weight;
      }
      break;
    }
  }

  TLorentzVector nuP4 = *(decayer->GetDecay(0));
  if (nthrows >= nmaxthrows) {

    *Decay_weight = weight / wmax;
    nuP4.Boost(parentP4.BoostVector());

  }
  else {
    *Decay_weight = 1.;

    nuP4.Boost(parentP4.BoostVector());
  }

  return nuP4;
}

//___________________________________________________________________________
std::vector<double> Calc_HSN_rest_angs(double ang_lab, double beta_meson, double beta_HSN_rest) {
  double tan_theta = TMath::Tan(ang_lab);

  double gamma_meson = (1. / (TMath::Sqrt(1 - (beta_meson * beta_meson))));


  double beta_ratio = beta_meson / beta_HSN_rest;
  double a = -1 * (tan_theta * gamma_meson);
  double b = 1;
  double c = beta_ratio * (tan_theta * gamma_meson);


  double t1 = TMath::Sqrt(pow(a, 2) + pow(b, 2) - pow(c, 2));

  double s1 = (b + t1) / (a + c);
  double s2 = (b - t1) / (a + c);


  std::vector<double> angs_rest;
  angs_rest.clear();
  angs_rest.push_back(2 * TMath::ATan(s1));
  angs_rest.push_back(2 * TMath::ATan(s2));

  return angs_rest;

}

//___________________________________________________________________________
std::vector<double> HSNangfactor(double ang_lab, double beta_meson, double beta_HSN_rest) {

  std::vector<double> rest_angs;
  std::vector<double> angfacts;
  rest_angs.clear();
  angfacts.clear();

  double ang_lab_rad = TMath::ACos(ang_lab);
  if (debug1) std::cout << "ang_lab_rad = " << ang_lab_rad << std::endl;
  rest_angs = Calc_HSN_rest_angs( TMath::ACos(ang_lab), beta_meson, beta_HSN_rest);
  if (debug1) std::cout << "rest_angs = " << rest_angs[0] << "," << rest_angs[1] << std::endl;
  double gamma_meson = (1 / (TMath::Sqrt(1 - (beta_meson * beta_meson))));
  if (debug1) std::cout << "gamma_meson = " << gamma_meson << std::endl;
  double beta_ratio = beta_meson / beta_HSN_rest;
  if (debug1) std::cout << "beta_ratio = " << beta_ratio << std::endl;

  double t1 = 1 - (TMath::Power(TMath::Tan(ang_lab_rad), 2) * TMath::Power(gamma_meson, 2) * (TMath::Power(beta_ratio, 2) - 1));
  if (debug1) std::cout << "t1 = " << t1 << std::endl;
  double t2_1 = gamma_meson * (1 - (beta_ratio / TMath::Sqrt(t1)));
  if (debug1) std::cout << "t2_1 = " << t2_1 << std::endl;
  double t2_2 = gamma_meson * (1 + (beta_ratio / TMath::Sqrt(t1)));
  if (debug1) std::cout << "t2_2 = " << t2_2 << std::endl;
  double t3 = 1 + (TMath::Power(TMath::Tan(ang_lab_rad), 2) * TMath::Power(gamma_meson, 2));
  if (debug1) std::cout << "t3 = " << t3 << std::endl;


  double w1 = TMath::Abs((t2_1 / t3) * (TMath::Sin(rest_angs[0]) / TMath::Sin(ang_lab_rad)) * (1 / TMath::Power(ang_lab, 2)));
  double w2 = TMath::Abs((t2_2 / t3) * (TMath::Sin(rest_angs[1]) / TMath::Sin(ang_lab_rad)) * (1 / TMath::Power(ang_lab, 2)));

  angfacts.push_back(w1);
  angfacts.push_back(w2);

  return angfacts;

}


//___________________________________________________________________________
int calcHSNWgt(bsim::Dk2Nu* dk2nu, const TVector3& xyz, std::vector<double>& eHSN, std::vector<double>& HSN_wgts_xy, double HSN_mass, int Ndecay_mod) {

//   const double kPIMASS = 0.13957;
//   const double kKMASS  = 0.49368;
//   const double kK0MASS = 0.49767;
//   const double kMUMASS = 0.105658389;
//   const double kOMEGAMASS = 1.67245;

//   const int kpdg_nue       =   12;  // extended Geant 53
//   const int kpdg_nuebar    =  -12;  // extended Geant 52
//   const int kpdg_eplus     =    -11;  // Geant  5
//   const int kpdg_eminus    =     11;  // Geant  6
//   const int kpdg_muplus     =   -13;  // Geant  5
//   const int kpdg_muminus    =    13;  // Geant  6
//   const int kpdg_pionplus   =   211;  // Geant  8
//   const int kpdg_pionminus  =  -211;  // Geant  9
//   const int kpdg_pi0        =   111;  // Geant  9
//   const int kpdg_k0long     =   130;  // Geant 10  ( K0=311, K0S=310 )
// // const int kpdg_k0short    =   310;  // Geant 16 //not used none of the dk2nu files
// // const int kpdg_k0mix      =   311;  //not used
//   const int kpdg_kaonplus   =   321;  // Geant 11
//   const int kpdg_kaonminus  =  -321;  // Geant 12
// // const int kpdg_omegaminus =  3334;  // Geant 24 dont get used

//   const int kpdg_muplus_G     =   5;  // Geant  5
//   const int kpdg_muminus_G    =    6;  // Geant  6
//   const int kpdg_pionplus_G   =   8;  // Geant  8
//   const int kpdg_pionminus_G  =  9;  // Geant  9
//   const int kpdg_k0long_G     =   10;  // Geant 10  ( K0=311, K0S=310 )
//   const int kpdg_k0short_G    =   16;  // Geant 16
//   const int kpdg_k0mix_G      =   311;
//   const int kpdg_kaonplus_G   =   11;  // Geant 11
//   const int kpdg_kaonminus_G  =  12;  // Geant 12
//   const int kpdg_omegaminus_G =  24;  // Geant 24
//   const int kpdg_omegaplus_G  = 32;  // Geant 32

  const double kRDET = 100.0;   // set to flux per 100 cm radius

  double xpos = xyz.X();
  double ypos = xyz.Y();
  double zpos = xyz.Z();




  // in principle we should get these from the particle DB
  //  but for consistency testing use the hardcoded values

  double parent_mass = kPIMASS;
  switch ( dk2nu->decay.ptype ) {
  case kpdg_pionplus:
  case kpdg_pionminus:
    parent_mass = kPIMASS;
    break;
  case kpdg_kaonplus:
  case kpdg_kaonminus:
    parent_mass = kKMASS;
    break;
  case kpdg_k0long:
  case kpdg_k0short:
  case kpdg_k0mix:
    parent_mass = kK0MASS;
    break;
  case kpdg_muplus:
  case kpdg_muminus:
    parent_mass = kMUMASS;
    break;
  case kpdg_omegaminus:
  case kpdg_omegaplus:
    parent_mass = kOMEGAMASS;
    break;
  default:
    std::cerr << "bsim::calcEnuWgt unknown particle type " << dk2nu->decay.ptype
              << std::endl << std::flush;
    //assert(0);
    return 1;
  }

  if (debug1) std::cout << "decay->pdpz: " << dk2nu->decay.pdpz << std::endl;
  double parentp2 = ( dk2nu->decay.pdpx * dk2nu->decay.pdpx +
                      dk2nu->decay.pdpy * dk2nu->decay.pdpy +
                      dk2nu->decay.pdpz * dk2nu->decay.pdpz );
  // double parent_energy = TMath::Sqrt( parentp2 +
  //                                    parent_mass*parent_mass);
  double parentp = TMath::Sqrt( parentp2 );


  // Get angle from parent line of flight to chosen point in beam frame
  double rad = TMath::Sqrt( (xpos - dk2nu->decay.vx) * (xpos - dk2nu->decay.vx) +
                            (ypos - dk2nu->decay.vy) * (ypos - dk2nu->decay.vy) +
                            (zpos - dk2nu->decay.vz) * (zpos - dk2nu->decay.vz) );

  if (debug1) std::cout << "rad from parent line to point =  " << rad << std::endl;



  const std::map<int, int> mother_pdgs = {
    { 1,  kpdg_k0long    },
    { 2,  kpdg_k0long    },
    { 3,  kpdg_k0long    },
    { 4,  kpdg_k0long    },
    { 5,  kpdg_kaonplus  },
    { 6,  kpdg_kaonplus  },
    { 7,  kpdg_kaonplus  },
    { 8,  kpdg_kaonminus },
    { 9,  kpdg_kaonminus },
    { 10, kpdg_kaonminus },
    { 11, kpdg_muplus    },
    { 12, kpdg_muminus   },
    { 13, kpdg_pionplus  },
    { 14, kpdg_pionminus },
    { 15, kpdg_kaonplus },
    { 16, kpdg_kaonminus },
    { 17, kpdg_pionplus },
    { 18, kpdg_pionminus }
  };

  const std::map<int, std::vector<int>> daughter_pdgs = {
    { 1,  { kpdg_pionminus, kpdg_eplus    } },
    { 2,  { kpdg_pionplus,  kpdg_eminus   } },
    { 3,  { kpdg_pionminus, kpdg_muplus   } },
    { 4,  { kpdg_pionplus,  kpdg_muminus  } },
    { 5,  { kpdg_muplus                   } },
    { 6,  { kpdg_pi0,       kpdg_eplus    } },
    { 7,  { kpdg_pi0,       kpdg_muplus   } },
    { 8,  { kpdg_muminus                  } },
    { 9,  { kpdg_pi0,       kpdg_eminus   } },
    { 10, { kpdg_pi0,       kpdg_muminus  } },
    { 11, { kpdg_nue,       kpdg_eplus    } },
    { 12, { kpdg_nuebar,    kpdg_eminus   } },
    { 13, { kpdg_muplus                   } },
    { 14, { kpdg_muminus                  } },
    { 15, { kpdg_eplus                    } },
    { 16, { kpdg_eminus                   } },
    { 17, { kpdg_eplus                    } },
    { 18, { kpdg_eminus                   } }
  };

// based on Geant4: G4KL3DecayChannel.cc

  const std::map<int, double> pLambda_map = {
    { 6, 0.0286 }, { 9,  0.0286 }, // k+- -> e+-
    { 7, 0.033  }, { 10, 0.033  }, // k+- -> mu+-
    { 1, 0.0300 }, { 2,  0.0300 }, // K0L -> e+-
    { 3, 0.034  }, { 4,  0.034  }  // K0L -> mu+-
  };
  const std::map<int, double> pXi0_map = {
    { 6, -0.35 }, { 9,  -0.35 }, // k+- -> e+-
    { 7, -0.35 }, { 10, -0.35 }, // k+- -> mu+-
    { 1, -0.11 }, { 2,  -0.11 }, // K0L -> e+-
    { 3, -0.11 }, { 4,  -0.11 }  // K0L -> mu+-
  };

  //Owen's stuff starts
  const std::map<int, std::vector<int>>& dpdgs = daughter_pdgs;
  const std::map<int, double> Qv = get_Q_vals(dpdgs, mother_pdgs);


  //const int decay_type = dk2nu->decay.ndecay;


  //parent mass already defined
  //cout<<decay_type<<endl;
  if (Qv.find(Ndecay_mod) == Qv.end()) {
    std::cout << "Unknown decay: " << Ndecay_mod << std::endl;
    return 4;
  }


  // end of if on decay not possible for this HSN_mass return
  if (Qv.at(Ndecay_mod) < HSN_mass) {
    if (debug1) std::cout << "decay not poss Qvalue" << std::endl;
    return 2;
  }


  //decay is possible

  double Decay_weight = 0;
  TLorentzVector parentP4;
  parentP4.SetXYZM(dk2nu->decay.pdpx, dk2nu->decay.pdpy, dk2nu->decay.pdpz, parent_mass);
  TLorentzVector hnu = make_decay(pLambda_map, pXi0_map, get_decayer(Ndecay_mod, HSN_mass, dpdgs, mother_pdgs), parentP4,  Ndecay_mod, &Decay_weight);
  // Get solid angle/4pi for detector element
  //     // small angle approximation, fixed by Alex Radovic
  //         //SAA//  double sangdet = ( kRDET*kRDET /
  //             //SAA//                   ( (zpos-decay->Vz)*(zpos-decay->Vz) ) ) / 4.0;

  hnu.Boost(-parentP4.BoostVector());
  double sanddetcomp = TMath::Sqrt( ( (xpos - dk2nu->decay.vx) * (xpos - dk2nu->decay.vx) ) +
                                    ( (ypos - dk2nu->decay.vy) * (ypos - dk2nu->decay.vy) ) +
                                    ( (zpos - dk2nu->decay.vz) * (zpos - dk2nu->decay.vz) )   );

  double sangdet = (1.0 - TMath::Cos(TMath::ATan( kRDET / sanddetcomp ))) / 2.0;


  const TVector3 pos(xpos - dk2nu->decay.vx, ypos - dk2nu->decay.vy, zpos - dk2nu->decay.vz);

  if ( parentp == 0 ) { //at rest
    eHSN.push_back(hnu.E());
    HSN_wgts_xy.push_back(Decay_weight * sangdet);
    // if (debug1) cout << hnu.E() << endl;
    return 5;
  }


  //not at rest

  double Betaratio = parentP4.Beta() / hnu.Beta();



  double costh_pardet = -999., theta_pardet = -999.;
  costh_pardet = ( dk2nu->decay.pdpx * (xpos - dk2nu->decay.vx) +
                   dk2nu->decay.pdpy * (ypos - dk2nu->decay.vy) +
                   dk2nu->decay.pdpz * (zpos - dk2nu->decay.vz) )
                 / ( parentp * rad);


  // if ( costh_pardet >  1.0 ) costh_pardet =  1.0;   \\Weighting functions blows up for cos=1 but tends to a nice value as you approach it
  // if ( costh_pardet < -1.0 ) costh_pardet = -1.0;

  if ( costh_pardet >  1.0 ) costh_pardet =  0.999999999999;
  if ( costh_pardet < -1.0 ) costh_pardet = -0.999999999999;

  theta_pardet = TMath::ACos(costh_pardet);

  std::vector<double> hsnangfacts = HSNangfactor(costh_pardet, parentP4.Beta(), hnu.Beta());





  TLorentzVector nhnm1;
  TLorentzVector nhnm2;
  TVector3 cm_dir = (1. / parentP4.Vect().Mag()) * parentP4.Vect();
  std::vector<double> rest_angs = Calc_HSN_rest_angs(theta_pardet, parentP4.Beta(), hnu.Beta());
  nhnm1.SetXYZM(hnu.Vect().Mag()*TMath::Sin(rest_angs[0]), 0, hnu.Vect().Mag()*TMath::Cos(rest_angs[0]), hnu.M());
  nhnm1.RotateUz(cm_dir);
  nhnm2.SetXYZM(hnu.Vect().Mag()*TMath::Sin(rest_angs[1]), 0, hnu.Vect().Mag()*TMath::Cos(rest_angs[1]), hnu.M());
  nhnm2.RotateUz(cm_dir);


  nhnm1.Boost(parentP4.BoostVector());
  nhnm2.Boost(parentP4.BoostVector());


  if (parentP4.Beta() / hnu.Beta() > 1) { //two solutions

    double maxang = TMath::ATan(TMath::Sqrt(1. / (parentP4.Gamma() * parentP4.Gamma() * ((Betaratio * Betaratio) - 1)) )); //calcs the max acheivable angle of HNL from parent direction
    if ((theta_pardet) > maxang) { //max acheivable angle, angle func will have (correctly) failed
      if (debug1) std::cout << "NOT possible to reach detector" << std::endl;
      return 1;

    }

    if ((hsnangfacts[0] > 100000) | (hsnangfacts[1] > 100000)) { //pciked by eye
      std::cout << "Found Large Weight, occurs when angle to TPC approaches max ang " << std::endl;
      std::cout << "Large weights: " << hsnangfacts[0] << "," << hsnangfacts[1] << std::endl;
      std::cout << "theta_pardet:" << theta_pardet << std::endl;
      std::cout << "max_ang:" << maxang << std::endl;
      return 8;
    }

    HSN_wgts_xy.push_back(Decay_weight * sangdet *  hsnangfacts[0] );
    HSN_wgts_xy.push_back(Decay_weight * sangdet *  hsnangfacts[1] );
    eHSN.push_back(nhnm1.E());
    eHSN.push_back(nhnm2.E());


  }//end two sols

  else if (parentP4.Beta() / hnu.Beta() < 1) { //only one solution in theta range 0 to pi

    if (theta_pardet < TMath::Pi() / 2) {
      eHSN.push_back(nhnm2.E());
      HSN_wgts_xy.push_back(Decay_weight * sangdet *  hsnangfacts[1] );
      if (rest_angs[0] > 0) { //if the ang and wight function are correct this should never trigger
        std::cout << "seem to be chucking away a good sol, make sure ang func and weigh fucntion match up" << std::endl;
        exit (EXIT_FAILURE);
      }

    }

    if (theta_pardet > TMath::Pi() / 2) {
      eHSN.push_back(nhnm1.E());
      HSN_wgts_xy.push_back(Decay_weight * sangdet *  hsnangfacts[0] );
      if (rest_angs[1] > 0) { //if the ang and wight function are correct this should never trigger
        std::cout << "seem to be chucking away a good sol, make sure ang func and weigh fucntion match up" << std::endl;
        exit (EXIT_FAILURE);
      }
    }

  } //end one sols

  else {
    std::cout << "weird beta_ratio" << std::endl;
    exit (EXIT_FAILURE);
  }




  if (Ndecay_mod == 11 || Ndecay_mod == 12) {
    // TODO: muon decay anisotropy
  }

  return 0;
}



// #endif

// //#ifdef BNB_HNLFlux_funcs_cxx

// //___________________________________________________________________________
// BNB_HNLFlux::BNB_HNLFlux(string pattern) {

//   const char* path = "/uboone/app/users/ogoodwin/BNB_HNLFlux";
//   if ( path ) {
//     TString libs = gSystem->GetDynamicPath();
//     libs += ":";
//     libs += path;
//     //libs += "/lib";
//     gSystem->SetDynamicPath(libs.Data());
//     gSystem->Load("FluxNtuple_C.so");
//   }

//   cflux = new TChain("h10");
//   cflux->Add(pattern.c_str());

//   Nfiles = cflux->GetNtrees();
//   cout << "Number of files: " << Nfiles << endl;

//   //Inizialise histos
//   TString titleBase1 = "Neutrino Flux;";
//   TString titleBase2 = " Energy [GeV];";
//   TString titleBase3 = " / cm^{2} / 6e20 POT";


//   numuFluxHisto = new TH1D("numuFluxHisto", (titleBase1 + "#nu_{#mu}" + titleBase2 + "#nu_{#mu}" + titleBase3), histNbins, histMin, histMax);
//   anumuFluxHisto = new TH1D("anumuFluxHisto", (titleBase1 + "#bar{#nu}_{#mu}" + titleBase2 + "#bar{#nu}_{#mu}" + titleBase3), histNbins, histMin, histMax);
//   nueFluxHisto = new TH1D("nueFluxHisto", (titleBase1 + "#nu_{e}" + titleBase2 + "#nu_{e}" + titleBase3), histNbins, histMin, histMax);
//   anueFluxHisto = new TH1D("anueFluxHisto", (titleBase1 + "#bar{#nu}_{e}" + titleBase2 + "#bar{#nu}_{e}" + titleBase3), histNbins, histMin, histMax);



//   numuCCHisto = new TH1D("numuCCHisto", "numu CC; #nu_{#mu} Energy [GeV]; #nu_{#mu} CC / 79 ton / 6e20 POT", histNbins, histMin, histMax);

//   HSNmuEvDist = new TH2D("HSNmuEvDist", "distHist; N_{#mu} Energy [GeV];Dist [cm]", histNbins, histMin, histMax, DistNbins, DistMin, DistMax);
//   numuEvDist = new TH2D("nummuEvDist", "distHist; #nu_{#mu}  Energy [GeV];Dist [cm]", histNbins, histMin, histMax, DistNbins, DistMin, DistMax);
//   numuEvDistCC = new TH2D("nummuEvDistCC", "distHist; #nu_{#mu} Energy [GeV];Dist [cm]", histNbins, histMin, histMax, DistNbins, DistMin, DistMax);
//   Decay_types = new TH1D("Decay_types", "Decay_types; Decay_types; #nu_{#mu} CC / 79 ton / 6e20 POT", 14, 0.5, 14.5);


// }

// BNB_HNLFlux::~BNB_HNLFlux() {

// }


//#endif // #ifdef BNB_HNLFlux_cxx
