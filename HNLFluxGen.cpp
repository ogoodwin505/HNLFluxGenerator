#include "TROOT.h"
#include "TRandom.h"
#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH3F.h"
#include "TProfile2D.h"
#include "TChain.h"
#include "TStyle.h"
#include "TString.h"
#include <vector>
#include <list>
#include <algorithm>
#include <iostream>
#include "TMath.h"

#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"
// #include "Parameter.cpp"
#include "Parameter.h"
#include "HNLFlux_funcs.h"
#include "HNLFluxGen.h"




static const int numu  = 14;
static const int anumu = -14;
static const int nue   = 12;
static const int anue  = -12;
bool debug = false;

int main(int argc, char* argv[]){

  double HSN_mass = 0.260;
  
  string parameter_file;
  if(argc==2)
    parameter_file = argv[1];
  else
    parameter_file = "MicroBooNEConfig.dat";
  std::cout<<"Using Par File: "<<parameter_file<<std::endl;
  Parameter par(parameter_file);
  std::string outputloc=par.Output_loc();
  std::cout<<outputloc<<std::endl;
  TString filename;
  filename.Form("HNLFlux_mass%.3fGeV.root", HSN_mass);

  TString filepath=outputloc+filename;
  //65921088
  std::cout << "Simulating Flux for HNL with mass: " << HSN_mass << "GeV" << std::endl;
  std::cout << filename << std::endl;
  TFile* f = new TFile(filepath, "RECREATE");


  TChain *cflux = new TChain("dk2nuTree");

  cflux->Add("dk2nuFiles/nubeam*.dk2nu.root");


  bsim::Dk2Nu*  dk2nu  = new bsim::Dk2Nu;
  cflux->SetBranchAddress("dk2nu",&dk2nu);


  TString titleBase1 = "Neutrino Flux;";
  TString titleBase2 = " Energy [GeV] ;";
  TString titleBase3 = " / cm^{2} / POT";
  int Ntype_mod;
  int Ndecay_mod;
  double ratio_K_plus_mu_e;
  double ratio_K_minus_mu_e;
  double ratio_pi_plus_mu_e;
  double ratio_pi_minus_mu_e;
  double scale_K_plus_e;
  double scale_K_minus_e;
  double scale_pi_plus_e;
  double scale_pi_minus_e;

  for (int decay_type = 0; decay_type < 19; ++decay_type) {
    numu_decay_type_hists[decay_type] = new TH1D( Form("numu_Flux_decaytype_%d", decay_type), Form("numu_Flux_decaytype_%d", decay_type) , histNbins, histMin, histMax);
    anumu_decay_type_hists[decay_type] = new TH1D( Form("anumu_Flux_decaytype_%d", decay_type), Form("anumu_Flux_decaytype_%d", decay_type) , histNbins, histMin, histMax);
    nue_decay_type_hists[decay_type] = new TH1D( Form("nue_Flux_decaytype_%d", decay_type), Form("nue_Flux_decaytype_%d", decay_type) , histNbins, histMin, histMax);
    anue_decay_type_hists[decay_type] = new TH1D( Form("anue_Flux_decaytype_%d", decay_type), Form("anue_Flux_decaytype_%d", decay_type) , histNbins, histMin, histMax);

    HSNmu_decay_type_hists[decay_type] = new TH1D( Form("HSNmu_Flux_decaytype_%d", decay_type), Form("HSNmu_Flux_decaytype_%d", decay_type) , histNbins, histMin, histMax);
    aHSNmu_decay_type_hists[decay_type] = new TH1D( Form("aHSNmu_Flux_decaytype_%d", decay_type), Form("aHSNmu_Flux_decaytype_%d", decay_type) , histNbins, histMin, histMax);
    HSNe_decay_type_hists[decay_type] = new TH1D( Form("HSNe_Flux_decaytype_%d", decay_type), Form("HSNe_Flux_decaytype_%d", decay_type) , histNbins, histMin, histMax);
    aHSNe_decay_type_hists[decay_type] = new TH1D( Form("aHSNe_Flux_decaytype_%d", decay_type), Form("aHSNe_Flux_decaytype_%d", decay_type) , histNbins, histMin, histMax);

  }

  // HSNmu
  HSNmuFluxHisto = new TH1D("HSNmuFluxHisto", (titleBase1 + "N_{#mu}" + titleBase2 + "N_{#mu}" + titleBase3), histNbins, histMin, histMax);
  // aHSNmu
  aHSNmuFluxHisto = new TH1D("aHSNmuFluxHisto", (titleBase1 + "#bar{N}_{#mu}" + titleBase2 + "#bar{N}_{#mu}" + titleBase3), histNbins, histMin, histMax);
  // HSNe
  HSNeFluxHisto = new TH1D("HSNeFluxHisto", (titleBase1 + "N_{e}" + titleBase2 + "N_{e}" + titleBase3), histNbins, histMin, histMax);
  // aHSNe
  aHSNeFluxHisto = new TH1D("aHSNeFluxHisto", (titleBase1 + "#bar{N}_{e}" + titleBase2 + "#bar{N}_{e}" + titleBase3), histNbins, histMin, histMax);

  WeightsHisto = new TH1D("WeightsHisto", (titleBase1 + "Position weight"), 500, 0.0, 1.0);

  numuFluxHisto = new TH1D("numuFluxHisto", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
  anumuFluxHisto = new TH1D("anumuFluxHisto", (titleBase1 + "#bar{#nu}_{#mu}" + titleBase2 +"#bar{#nu}_{#mu}" + titleBase3),histNbins,histMin,histMax);
  nueFluxHisto = new TH1D("nueFluxHisto", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
  anueFluxHisto = new TH1D("anueFluxHisto", (titleBase1 + "#bar{#nu}_{e}" + titleBase2 + "#bar{#nu}_{e}" + titleBase3),histNbins,histMin,histMax);
  numuCCHisto = new TH1D("numuCCHisto", "numu CC; #nu_{#mu} Energy [GeV]; #nu_{#mu} CC / 79 ton / 6e20 POT", histNbins, histMin, histMax);
  HSNmuEvDist = new TH2D("HSNmuEvDist", "distHist; N_{#mu} Energy [GeV];Dist [cm]", histNbins, histMin, histMax, DistNbins, DistMin, DistMax);
  numuEvDist = new TH2D("nummuEvDist", "distHist; #nu_{#mu}  Energy [GeV];Dist [cm]", histNbins, histMin, histMax, DistNbins, DistMin, DistMax);
  numuEvDistCC = new TH2D("nummuEvDistCC", "distHist; #nu_{#mu} Energy [GeV];Dist [cm]", histNbins, histMin, histMax, DistNbins, DistMin, DistMax);
  Decay_types = new TH1D("Decay_types", "Decay_types; Decay_types; #nu_{#mu} CC / 79 ton / 6e20 POT", 14, 0.5, 14.5);


  Long64_t nflux = cflux->GetEntries();
  std::cout << "Total number of entries: " << nflux << std::endl;
  int totmeson = 0;
  int Qval = 0;
  int reachDec = 0;
  int missDec = 0;
  int DAR = 0;
  int nans = 0;
  int large_w = 0;
  AccumulatedPOT = 0;

  int rand_e_mu = 0;

  //Integers to store the number of different K and pi decays
  int K_like_count = 0;
  int pi_like_count = 0;
  int K_plus_mu_count = 0;
  int K_plus_e_count = 0;
  int K_minus_mu_count = 0;
  int K_minus_e_count = 0;
  int pi_plus_mu_count = 0;
  int pi_plus_e_count = 0;
  int pi_minus_mu_count = 0;
  int pi_minus_e_count = 0;
  int K_plus_tot;
  int K_minus_tot;
  int pi_plus_tot;
  int pi_minus_tot;



  for (Long64_t i=0; i < nflux; ++i ) {
    cflux->GetEntry(i);
    Ntype_mod = dk2nu->decay.ntype;
    Ndecay_mod = dk2nu->decay.ndecay;
    // Inlcuding the electron like decays

    //AccumulatedPOT += dk2nu->potnum;

    //std::cout << "iteration is " << i << std::endl;

    rand_e_mu = Random_e_mu(i);

    switch (dk2nu->decay.ndecay) {
    case 5:
      K_like_count++;
      if (rand_e_mu == 2) {
        Ntype_mod = numu;
        K_plus_mu_count++;
      }
      else {
        Ntype_mod = nue;
        K_plus_e_count++;
        Ndecay_mod = 15;
      }
      break;
    case 13:
      pi_like_count++;
      if (rand_e_mu == 2) {
        Ntype_mod = numu;
        pi_plus_mu_count++;
      }
      else {
        Ntype_mod = nue;
        pi_plus_e_count++;
        Ndecay_mod = 17;
      }
      break;
    case 8:
      K_like_count++;
      if (rand_e_mu == 2) {
        Ntype_mod = anumu;
        K_minus_mu_count++;
      }
      else {
        Ntype_mod = anue;
        K_minus_e_count++;
        Ndecay_mod = 16;
      }
      break;
 
    case 14:
      pi_like_count++;
      if (rand_e_mu == 2) {
        Ntype_mod = anumu;
        pi_minus_mu_count++;
      }
      else {
        Ntype_mod = anue;
        pi_minus_e_count++;
        Ndecay_mod = 18;
      }
      break;

      case 999: //invalid decay 
        continue;
    }

    // Alert the user
    if (i % 100000 == 0) std::cout << "On entry " << i << std::endl;
    if (treeNumber != cflux->GetTreeNumber()) {
      treeNumber = cflux->GetTreeNumber();
      std::cout << "Moving to tree number " << treeNumber << "." << std::endl;
      // AccumulatedPOT += estimate_pots(highest_evtno);
      AccumulatedPOT += 100000; //per file
      std::cout << "AccumulatedPOT: " << AccumulatedPOT << std::endl;
    }

    std::vector<double> HSN_wgts_xy = {}; // HSN weights
    std::vector<double> eHSN = {}; // HSN energies
    double wgt_xy = 0.;  // neutrino weight
    double enu    = 0.;  // neutrino energy in lab frame

    // Pick a random point in the TPC (in detector coordinates)
    TVector3 xyz_det = RandomInTPC(i,par);
    if (debug) std::cout << "xyz_det = [" << xyz_det.X() << ", " << xyz_det.Y() << ", " << xyz_det.Z() << "]" << std::endl;

    // From detector to beam coordinates
    TVector3 xyz_beam = FromDetToBeam(xyz_det,par);
    if (debug) std::cout << "xyz_beam = [" << xyz_beam.X() << ", " << xyz_beam.Y() << ", " << xyz_beam.Z() << "]" << std::endl;

    // Calculate the weight for neutrino
    int ret = calcEnuWgt(dk2nu, xyz_beam, enu, wgt_xy);

    // Calculate the weight for HSN
    if (debug) std::cout << "running here?" << std::endl;
    int ret1 = calcHSNWgt(dk2nu, xyz_beam, eHSN, HSN_wgts_xy, HSN_mass, Ndecay_mod);
    if (debug) std::cout << "ret1: " << ret1 << std::endl;
    //Find distance between meson decay vertex and detector
    double Pathlen = findPathLen(dk2nu, xyz_beam);
    if (debug) std::cout << "ndecay" << dk2nu->decay.ndecay << std::endl;
    if (debug) std::cout << "ntype" << dk2nu->decay.ntype << std::endl;
    if (debug) std::cout << "nimwt=" << dk2nu->decay.nimpwt << std::endl;
    totmeson++;
    if (ret1 == 0) reachDec++;
    if (ret1 == 1) missDec++;
    if (ret1 == 2) Qval++;
    if (ret1 == 5) DAR = DAR + dk2nu->decay.nimpwt;
    if (ret1 ==8)  large_w++;
    if (ret != 0) std::cout << "Error with calcEnuWgt. Return " << ret << std::endl;

    if (debug) std::cout << "wgt_xy " << wgt_xy << std::endl;

    double weight = wgt_xy * dk2nu->decay.nimpwt * fDefaultWeightCorrection;
    for (unsigned int j = 0; j < HSN_wgts_xy.size(); ++j) {
      WeightsHisto->Fill(HSN_wgts_xy[j]);
      // if (HSN_wgts_xy[j] > 0.10) {
      //   std::cout << "found big w in it: " << std::endl;
      //   std::cout << i << std::endl;
      //   std::cout << "Large weight: " << HSN_wgts_xy[j] << std::endl;
      //   ++;
      //   HSN_wgts_xy[j] = 0;
      //   exit (EXIT_FAILURE);

      // }
    }
    for (unsigned int j = 0; j < HSN_wgts_xy.size(); ++j) { //does this not run if unedited
      if (std::isnan(HSN_wgts_xy[j])) {
        std::cout << "found nan in weight in it: " << std::endl;
        std::cout << i << std::endl;
        nans++;
        HSN_wgts_xy[j] = 0;
        exit (EXIT_FAILURE);
      }
      Decay_types->Fill(Ndecay_mod, HSN_wgts_xy[j]* dk2nu->decay.nimpwt * fDefaultWeightCorrection);

    }
    for (unsigned int j = 0; j < eHSN.size(); ++j) { //does this not run if unedited
      if (std::isnan(eHSN[j])) {
        std::cout << "found nan energy in it: " << std::endl;
        std::cout << i << std::endl;
        HSN_wgts_xy[j] = 0;
        nans++;
        exit (EXIT_FAILURE);
      }
    }

    switch (Ntype_mod) {
    case numu:
      //std::cout<<enu<<","<<weight<<std::endl;
      numu_decay_type_hists[Ndecay_mod]->Fill(enu, weight);
      numuFluxHisto->Fill(enu, weight);
      numuEvDist->Fill(enu, Pathlen, weight);
      for (unsigned int j = 0; j < HSN_wgts_xy.size(); ++j) {

        HSNmu_decay_type_hists[Ndecay_mod]->Fill(eHSN[j], HSN_wgts_xy[j]* dk2nu->decay.nimpwt * fDefaultWeightCorrection);                                                  //does this not run if unedited
        HSNmuFluxHisto->Fill(eHSN[j], HSN_wgts_xy[j]* dk2nu->decay.nimpwt * fDefaultWeightCorrection);
        HSNmuEvDist->Fill(eHSN[j], Pathlen, HSN_wgts_xy[j]* dk2nu->decay.nimpwt * fDefaultWeightCorrection);


      }

      break;
    case anumu:

      anumu_decay_type_hists[Ndecay_mod]->Fill(enu, weight);
      anumuFluxHisto->Fill(enu, weight);

      for (unsigned int j = 0; j < HSN_wgts_xy.size(); ++j) { //does this not run if unedited
        aHSNmuFluxHisto->Fill(eHSN[j], HSN_wgts_xy[j]* dk2nu->decay.nimpwt * fDefaultWeightCorrection);
        aHSNmu_decay_type_hists[Ndecay_mod]->Fill(eHSN[j], HSN_wgts_xy[j]* dk2nu->decay.nimpwt * fDefaultWeightCorrection);
      }

      break;
    case nue:

      nue_decay_type_hists[Ndecay_mod]->Fill(enu, weight);
      nueFluxHisto->Fill(enu, weight);

      for (unsigned int j = 0; j < HSN_wgts_xy.size(); ++j) { //does this not run if unedited
        HSNeFluxHisto->Fill(eHSN[j], HSN_wgts_xy[j]* dk2nu->decay.nimpwt * fDefaultWeightCorrection);
        HSNe_decay_type_hists[Ndecay_mod]->Fill(eHSN[j], HSN_wgts_xy[j]* dk2nu->decay.nimpwt * fDefaultWeightCorrection);
      }

      break;
    case anue:
      anue_decay_type_hists[Ndecay_mod]->Fill(enu, weight);
      anueFluxHisto->Fill(enu, weight);

      for (unsigned int j = 0; j < HSN_wgts_xy.size(); ++j) { //does this not run if unedited
        aHSNeFluxHisto->Fill(eHSN[j], HSN_wgts_xy[j]* dk2nu->decay.nimpwt * fDefaultWeightCorrection);
        aHSNe_decay_type_hists[Ndecay_mod]->Fill(eHSN[j], HSN_wgts_xy[j]* dk2nu->decay.nimpwt * fDefaultWeightCorrection);
      }

      break;
    }

    // POT stuff
    if ( dk2nu->potnum > highest_evtno ) {
      highest_evtno = dk2nu->potnum;
    }
  } // end of loop over the entries

  //These store the total counts for different decays
  K_plus_tot = K_plus_mu_count + K_plus_e_count;
  K_minus_tot = K_minus_mu_count + K_minus_e_count;
  pi_plus_tot = pi_plus_mu_count + pi_plus_e_count;
  pi_minus_tot = pi_minus_mu_count + pi_minus_e_count;

  //These are used to scale the mu like histograms
  ratio_K_plus_mu_e = ((double)K_plus_tot) / K_plus_mu_count;
  ratio_K_minus_mu_e = ((double)K_minus_tot) / K_minus_mu_count;
  ratio_pi_plus_mu_e = ((double)pi_plus_tot) / pi_plus_mu_count;
  ratio_pi_minus_mu_e = ((double)pi_minus_tot) / pi_minus_mu_count;

  //These are used to scale the e like histograms
  scale_K_plus_e = (K_decay_e_fraction / K_decay_mu_fraction) * (((double)K_plus_tot) / K_plus_e_count);
  scale_K_minus_e = (K_decay_e_fraction / K_decay_mu_fraction) * (((double)K_minus_tot) / K_minus_e_count);
  scale_pi_plus_e = (pi_decay_e_fraction / pi_decay_mu_fraction) * (((double)pi_plus_tot) / pi_plus_e_count);
  scale_pi_minus_e = (pi_decay_e_fraction / pi_decay_mu_fraction) * (((double)pi_minus_tot) / pi_minus_e_count);
/*
  std::cout << "total Parents=" << totmeson<< "  With Non Allowed Qval = " << Qval << "  Can reach Decector (with non zero weight) = " << reachDec << std::endl;
  std::cout << "Can never reach Decector" << missDec  << "Decays at Rest = " << DAR << std::endl;
  std::cout<< "NaN weights or Energy (should be zero) = " << nans << std::endl;
  std::cout << "Number of Excluded 'large weights' = "<< large_w << std::endl;;
*/
  //***************************************
  //
  // POT scaling
  //
  //***************************************

  //AccumulatedPOT += estimate_pots(highest_evtno); // To account for last tree
  //AccumulatedPOT += 100000; //as last tree wont be counted
  double scale = NominalPOT / AccumulatedPOT;

  Decay_types -> Scale(scale);
  ///// Bin Width Scaleling, get units to per GeV

  float binwidth = (histMax - histMin) / histNbins; //width in GeV
  scale *= (1 / binwidth); //Scale to per 1GeV
  numuFluxHisto  -> Scale(scale);
  anumuFluxHisto -> Scale(scale);
  nueFluxHisto   -> Scale(scale);
  anueFluxHisto  -> Scale(scale);

  HSNmuFluxHisto  -> Scale(scale);
  aHSNmuFluxHisto -> Scale(scale);
  HSNeFluxHisto   -> Scale(scale);
  aHSNeFluxHisto  -> Scale(scale);
  numuEvDist -> Scale(scale);
  HSNmuEvDist -> Scale(scale);



  for (int decay_type = 0; decay_type < 19; ++decay_type) {

    numu_decay_type_hists[decay_type]->Scale(scale);
    anumu_decay_type_hists[decay_type]->Scale(scale);
    nue_decay_type_hists[decay_type]->Scale(scale);
    anue_decay_type_hists[decay_type]->Scale(scale);

    HSNmu_decay_type_hists[decay_type]->Scale(scale);
    aHSNmu_decay_type_hists[decay_type]->Scale(scale);
    HSNe_decay_type_hists[decay_type]->Scale(scale);
    aHSNe_decay_type_hists[decay_type]->Scale(scale);

    switch (decay_type) {
    case 5:
      numu_decay_type_hists[decay_type]->Scale(ratio_K_plus_mu_e);
      HSNmu_decay_type_hists[decay_type]->Scale(ratio_K_plus_mu_e);
      break;
    case 15:
      nue_decay_type_hists[decay_type]->Scale(scale_K_plus_e);
      HSNe_decay_type_hists[decay_type]->Scale(scale_K_plus_e);
      break;
    case 13:
      numu_decay_type_hists[decay_type]->Scale(ratio_pi_plus_mu_e);
      HSNmu_decay_type_hists[decay_type]->Scale(ratio_pi_plus_mu_e);
      break;
    case 17:
      nue_decay_type_hists[decay_type]->Scale(scale_pi_plus_e);
      HSNe_decay_type_hists[decay_type]->Scale(scale_pi_plus_e);
      break;
    case 8:
      anumu_decay_type_hists[decay_type]->Scale(ratio_K_minus_mu_e);
      aHSNmu_decay_type_hists[decay_type]->Scale(ratio_K_minus_mu_e);
      break;
    case 16:
      anue_decay_type_hists[decay_type]->Scale(scale_K_minus_e);
      aHSNe_decay_type_hists[decay_type]->Scale(scale_K_minus_e);
      break;
    case 14:
      anumu_decay_type_hists[decay_type]->Scale(ratio_pi_minus_mu_e);
      aHSNmu_decay_type_hists[decay_type]->Scale(ratio_pi_minus_mu_e);
      break;
    case 18:
      anue_decay_type_hists[decay_type]->Scale(scale_pi_minus_e);
      aHSNe_decay_type_hists[decay_type]->Scale(scale_pi_minus_e);
      break;
    }
  }
  std::cout << std::endl << ">>> TOTAL POT: " << AccumulatedPOT << std::endl << std::endl;


  //***************************************
  //
  // Apply now GENIE xsec
  //
  // source /nusoft/app/externals/setup
  // setup genie_xsec R-2_8_0   -q default
  // root -l  $GENIEXSECPATH/xsec_graphs.root
  // >  _file0->cd("nu_mu_Ar40")
  // >  tot_cc->Draw()
  //
  //***************************************

  const char* genieXsecPath = gSystem->ExpandPathName("$(GENIEXSECPATH)");
  if ( !genieXsecPath ) {
    std::cout << "$(GENIEXSECPATH) not defined." << std::endl;
    std::cout << "Please setup *genie_xsec*. (setup genie_xsec R-2_8_0   -q default)." << std::endl;
  }

  if ( genieXsecPath ) {
    TString genieXsecFileName = genieXsecPath;
    genieXsecFileName += "/xsec_graphs.root";
    TFile *genieXsecFile = new TFile(genieXsecFileName, "READ");
    genieXsecFile->cd("nu_mu_Ar40");
    genieXsecNumuCC = (TGraph *) gDirectory->Get("tot_cc");
    genieXsecFile->Close();

    // TSpline3* genieXsecSplineNumuCC = new TSpline3("genieXsecSplineNumuCC", genieXsecNumuCC, "", 0,6);

    double value;
    for (int i = 1; i < histNbins + 1; i++) {
      value = numuFluxHisto->GetBinContent(i);
      value *= genieXsecNumuCC->Eval(numuFluxHisto->GetBinCenter(i)); // Eval implies linear interpolation
      value *= (1e-38 * Ntarget / 40.); // 1/40 is due to I'm considering nu_mu_Ar40.

      numuCCHisto->SetBinContent(i, value);

      for (int j = 1; j < DistNbins + 1; j++) {
        value = numuEvDist->GetBinContent(i, j);
        value *= genieXsecNumuCC->Eval(numuEvDist->GetXaxis()->GetBinCenter(i)); // Eval implies linear interpolation
        value *= (1e-38 * Ntarget / 40.); // 1/40 is due to I'm considering nu_mu_Ar40.
        numuEvDistCC->SetBinContent(i, j, value);
      }

    }
  } // end if ( genieXsecPath )


  //***************************************
  //
  // Writing on file
  //
  //***************************************

  f->cd();


  for (int decay_type = 0; decay_type < 19; ++decay_type) {
    numu_decay_type_hists[decay_type]-> Write();
    anumu_decay_type_hists[decay_type]-> Write();
    nue_decay_type_hists[decay_type]-> Write();
    anue_decay_type_hists[decay_type]-> Write();

    HSNmu_decay_type_hists[decay_type]-> Write();
    aHSNmu_decay_type_hists[decay_type]-> Write();
    HSNe_decay_type_hists[decay_type]-> Write();
    aHSNe_decay_type_hists[decay_type]-> Write();
  }
  numuFluxHisto -> Draw("HIST");
  numuFluxHisto  -> Write();
  anumuFluxHisto -> Write();
  nueFluxHisto   -> Write();
  anueFluxHisto  -> Write();
  HSNmuFluxHisto  -> Write();
  aHSNmuFluxHisto -> Write();
  HSNeFluxHisto   -> Write();
  aHSNeFluxHisto  -> Write();
  HSNmuEvDist-> Write();
  numuEvDist-> Write();
  Decay_types->Write();

  WeightsHisto->Write();

  if ( genieXsecPath ) {
    numuCCHisto     -> Write();
    genieXsecNumuCC -> Write();
    numuEvDistCC -> Write();
  }

  f->Close();


  /*verx_hist->SetDirectory(fout);
  fout->Write();
  fout->Close();
 */
}
