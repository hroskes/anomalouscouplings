#include "ReducibleBackgroundAA_2015.C"

int didsetup = -1;
TFile *fFakeRates = 0;

bool ZXsetup(int production, TString dir) {
  if (didsetup == production) return true;
  didsetup = production;
  delete fFakeRates;
  if (production == 161221 || production == 170119) {
    fFakeRates = TFile::Open(dir+"/FakeRate_SS_2016D_hists.root");
    h1D_FRel_EB = (TH1F*)fFakeRates->Get("FR_SS_electron_EB");
    h1D_FRel_EE = (TH1F*)fFakeRates->Get("FR_SS_electron_EE");
    h1D_FRmu_EB = (TH1F*)fFakeRates->Get("FR_SS_muon_EB");
    h1D_FRmu_EE = (TH1F*)fFakeRates->Get("FR_SS_muon_EE");
  } else {
    return false;
  }
  return true;
}

