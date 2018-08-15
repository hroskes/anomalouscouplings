#include "ReducibleBackgroundAA_2015.C"

int didsetup = -1;
unique_ptr<TFile> fFakeRates = nullptr;

enum ZXsetupstatus {
  ZXsetupsuccess,
  ZXsetupbadproduction,
  ZXsetupfailed
};

ZXsetupstatus ZXsetup(int production, TString dir) {
  if (didsetup == production) return ZXsetupsuccess;
  didsetup = production;
  if (production == 170203 || production == 170222 || production == 170712 || production == 170825 || production == 180121 || production == 180224 || production == 180530 || production == 180721) {
    fFakeRates.reset(TFile::Open(dir+"/FakeRate_SS_Moriond368_hists.root"));
  } else if (production == 180416 || production == 180531 || production == 180722) {
    fFakeRates.reset(TFile::Open(dir+"/FakeRates_SS_Moriond18_hists.root"));
  } else {
    return ZXsetupbadproduction;
  }

  h1D_FRel_EB = dynamic_cast<TH1F*>(fFakeRates->Get("FR_SS_electron_EB"));
  h1D_FRel_EE = dynamic_cast<TH1F*>(fFakeRates->Get("FR_SS_electron_EE"));
  h1D_FRmu_EB = dynamic_cast<TH1F*>(fFakeRates->Get("FR_SS_muon_EB"));
  h1D_FRmu_EE = dynamic_cast<TH1F*>(fFakeRates->Get("FR_SS_muon_EE"));
  if (h1D_FRel_EB == 0 || h1D_FRel_EE == 0 || h1D_FRmu_EB == 0 || h1D_FRmu_EE == 0) {
    return ZXsetupfailed;
  }
  return ZXsetupsuccess;
}
