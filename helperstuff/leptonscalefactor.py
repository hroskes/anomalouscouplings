import ROOT

import CJLSTscripts
from utilities import cache, KeyDefaultDict

@cache
def LeptonSFHelper():
  return ROOT.LeptonSFHelper()

def fixleptonscalefactor(year, LepLepId, LepPt, LepEta, dataMCWeight):
  helper = LeptonSFHelper()
  updatedSF = (
    helper.getSF(year, LepLepId[0], LepPt[0], LepEta[0], LepEta[0], False) *
    helper.getSF(year, LepLepId[1], LepPt[1], LepEta[1], LepEta[1], False) *
    helper.getSF(year, LepLepId[2], LepPt[2], LepEta[2], LepEta[2], False) *
    helper.getSF(year, LepLepId[3], LepPt[3], LepEta[3], LepEta[3], False)
  )
  return updatedSF / dataMCWeight
