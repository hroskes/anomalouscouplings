import os

import ROOT

from utilities import KeyDefaultDict
from CJLSTscripts import CJLSTscriptsfolder, convertTGraphstoTH1Fs

class FakeRates(ROOT.FakeRates):
  def __init__(self, year):
    return super(FakeRates, self).__init__(
      {
         2016: os.path.join(CJLSTscriptsfolder, "FakeRate_SS_Moriond368.root"),
         2017: os.path.join(CJLSTscriptsfolder, "FakeRate_SS_Moriond18.root"),
         2018: os.path.join(CJLSTscriptsfolder, "FakeRate_SS_Moriond19.root"),
      }[year]
    )

__fakerates = KeyDefaultDict(FakeRates)

def getfakerate(year, leppt, lepeta, leplepid):
    return __fakerates[year].GetFakeRate(leppt, lepeta, leplepid)

from ROOT import CRZLLss, test_bit
