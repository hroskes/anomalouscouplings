#include "ReducibleBackgroundAA_2015.C"
#include "PyException.h"

int didsetup = -1;
TFile *fFakeRates = 0;

bool ZXsetup(int production, TString dir) {
  if (didsetup == production) return true;
  didsetup = production;
  delete fFakeRates;
  if (production == 170203) {
    fFakeRates = TFile::Open(dir+"/FakeRate_SS_Moriond368_hists.root");
    h1D_FRel_EB = dynamic_cast<TH1F*>(fFakeRates->Get("FR_SS_electron_EB"));
    h1D_FRel_EE = dynamic_cast<TH1F*>(fFakeRates->Get("FR_SS_electron_EE"));
    h1D_FRmu_EB = dynamic_cast<TH1F*>(fFakeRates->Get("FR_SS_muon_EB"));
    h1D_FRmu_EE = dynamic_cast<TH1F*>(fFakeRates->Get("FR_SS_muon_EEe"));
    if (h1D_FRel_EB == 0 || h1D_FRel_EE == 0 || h1D_FRmu_EB == 0 || h1D_FRmu_EE == 0) {
      throw PyException("Histograms not in the hists file or are not TH1Fs!");
    }
  } else {
    return false;
  }
  return true;
}

