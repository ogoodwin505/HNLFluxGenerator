#ifndef HNLFluxGen_h
#define HNLFluxGen_h

#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

#include "TChain.h"
#include "TSystem.h"
#include "TRandom.h"
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

#include "Parameter.h"

//#include "calcLocationWeights.h"


  //int Nfiles = 0;

  // static const int numu  = 56;
  // static const int anumu = 55;
  // static const int nue   = 53;
  // static const int anue  = 52;

  // static const int kpdg_nue_G       =   53;  // extended Geant 53
  // static const int kpdg_nuebar_G   =  52;  // extended Geant 52
  // static const int kpdg_numu_G      =   56;  // extended Geant 56
  // static const int kpdg_numubar_G   =  55;  // extended Geant 55

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
  
  double K_decay_mu_fraction = 0.6356;
  double K_decay_e_fraction = 0.00001582;
  double pi_decay_mu_fraction = 0.999877;
  double pi_decay_e_fraction = 0.000123;
 


  int highest_evtno = 0;
  double NominalPOT = 1;
  double fDefaultWeightCorrection = 1./(10000. * TMath::Pi());
  double Ntarget = 4.76e31/56.41e6*256.35*233*1036.8; //TPC active!!!
  double AccumulatedPOT=0.;
  int treeNumber = -1;

  double histMin = 0;
  double histMax = 5;
  int histNbins = 500;

  double DistMin=0;
  double DistMax=80000;
  int DistNbins= 50;



  //TChain *cflux;
  
  TH1D* numu_decay_type_hists[19];
  TH1D* anumu_decay_type_hists[19];
  TH1D* nue_decay_type_hists[19];
  TH1D* anue_decay_type_hists[19];



  TH1D* HSNmu_decay_type_hists[19];
  TH1D* aHSNmu_decay_type_hists[19];
  TH1D* HSNe_decay_type_hists[19];
  TH1D* aHSNe_decay_type_hists[19];
  //TH1D* decay_type_hists[14];
  //FluxNtuple *fluxNtuple;
  TH1D*  Decay_types;
  TH1D*  New_Decay_types;
  TH1D* numuFluxHisto;
  TH1D* numuFluxHisto_muminus;
  TH1D* numuFluxHisto_pionplus;
  TH1D* numuFluxHisto_kaonplus;
  TH1D* anumuFluxHisto;
  TH1D* nueFluxHisto;
  TH1D* anueFluxHisto;
  TH1D* numuCCHisto;

  TH1D* WeightsHisto;
  
  TGraph *genieXsecNumuCC;


  TH1D* HSNmuFluxHisto;
  TH1D* aHSNmuFluxHisto;
  TH1D* HSNeFluxHisto;
  TH1D* aHSNeFluxHisto;
  TH1D* HSNmuFluxHisto_kdar;
  TH1D* aHSNmuFluxHisto_kdar;
  TH2D* HSNmuEvDist;
  TH2D* numuEvDist;
  TH2D* numuEvDistCC;
  //const double HSN_mass=0.3;





  //BNB_HNLFlux(string pattern="/uboone/data/users/ogoodwin/BooNEG4BeamConv_tognumi/BNB_fluggfile_*.root");
//"/uboone/data/flux/numi/current/flugg_mn000z200i_20101117.gpcfgrid_lowth/flugg_mn000z200i_20101117.gpcfgrid_lowth_00*.root"
//
/*
  virtual ~BNB_HNLFlux();

  void CalculateFlux(double HSN_mass, std::string paramfile);
  TVector3 RandomInTPC(Long64_t counter, Parameter& par);
  TVector3 FromDetToBeam(const TVector3& det, Parameter& par);
  double estimate_pots(int highest_potnum);
  int calcEnuWgt( FluxNtuple* decay, const TVector3& xyz, double& enu, double& wgt_xy);

  map<int, double> get_Q_vals(const map<int, vector<int>>& daughters, const map<int, int> mother_pdgs);
  TGenPhaseSpace* get_decayer(int decay_type, double heavy_mass,const map<int, vector<int>> daughter_pdgs, const map<int, int> mother_pdgs);
  double DalitzDensity(const map<int, double> pLambda_map, const map<int, double> pXi0_map,const TLorentzVector& p4_k, const TLorentzVector& p4_pi, const TLorentzVector& p4_lep, const TLorentzVector& p4_Hnu, double decay_type);
  int calcHSNWgt( FluxNtuple* decay, const TVector3& xyz, vector<double>& eHSN, vector<double>& HSN_wgts_xy, double HSN_mass, int Ndecay_mod);
  TLorentzVector make_decay(const map<int, double> pLambda_map, const map<int, double> pXi0_map,TGenPhaseSpace* decayer, const TLorentzVector& parentP4,  int decay_type, double* Decay_weight);
  vector<double> Calc_HSN_rest_angs(double ang_lab,double beta_meson, double beta_HSN_rest);
  vector<double> HSNangfactor(double ang_lab, double beta_meson, double beta_HSN_rest);
  double findPathLen(FluxNtuple* decay,const TVector3& xyz);
  int Random_e_mu(Long64_t counter);
*/

#endif
