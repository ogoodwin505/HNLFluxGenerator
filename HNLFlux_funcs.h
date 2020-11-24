#ifndef HNLFlux_funcs_h
#define HNLFlux_funcs_h

#include <iostream>
#include <iomanip>
#include "TMath.h"

#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"

#include "Parameter.h"
#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"






//___________________________________________________________________________
int Random_e_mu(Long64_t counter);

// //___________________________________________________________________________
TVector3 RandomInTPC(Long64_t counter, Parameter& par);



//___________________________________________________________________________
TVector3 FromDetToBeam( const TVector3& det, Parameter& par);



//___________________________________________________________________________
double estimate_pots(int highest_potnum);



//___________________________________________________________________________
int calcEnuWgt(bsim::Dk2Nu* dk2nu, const TVector3& xyz,double& enu, double& wgt_xy);


//___________________________________________________________________________
double findPathLen(bsim::Dk2Nu* dk2nu, const TVector3& xyz);

//___________________________________________________________________________
double DalitzDensity(const std::map<int, double> pLambda_map, const std::map<int, double> pXi0_map, const TLorentzVector& p4_k, const TLorentzVector& p4_pi, const TLorentzVector& p4_lep, const TLorentzVector& p4_Hnu, double decay_type);

//___________________________________________________________________________
std::map<int, double> get_Q_vals(const std::map<int, std::vector<int>>& daughters, const std::map<int, int> mother_pdgs);


//___________________________________________________________________________
TGenPhaseSpace* get_decayer(int decay_type, double heavy_mass, const std::map<int, std::vector<int>> daughter_pdgs, const std::map<int, int> mother_pdgs);


//___________________________________________________________________________
TLorentzVector make_decay(const std::map<int, double> pLambda_map, const std::map<int, double> pXi0_map, TGenPhaseSpace* decayer, const TLorentzVector& parentP4,  int decay_type, double *Decay_weight);

//___________________________________________________________________________
std::vector<double> Calc_HSN_rest_angs(double ang_lab, double beta_meson, double beta_HSN_rest);

//___________________________________________________________________________
std::vector<double> HSNangfactor(double ang_lab, double beta_meson, double beta_HSN_rest);


//___________________________________________________________________________
int calcHSNWgt(bsim::Dk2Nu* dk2nu, const TVector3& xyz, std::vector<double>& eHSN, std::vector<double>& HSN_wgts_xy, double HSN_mass, int Ndecay_mod);



#endif

//#ifdef BNB_HNLFlux_funcs_cxx

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


//#endif 
