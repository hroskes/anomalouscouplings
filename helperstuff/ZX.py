import os

import ROOT

from utilities import KeyDefaultDict
from CJLSTscripts import CJLSTscriptsfolder, convertTGraphstoTH1Fs

CJLSTscriptsrelpath = os.path.relpath(CJLSTscriptsfolder)

class FakeRates(ROOT.FakeRates):
  def __init__(self, year):
    return super(FakeRates, self).__init__(
      {
         2016: os.path.join(CJLSTscriptsrelpath, "FakeRate_SS_Moriond368.root"),
         2017: os.path.join(CJLSTscriptsrelpath, "FakeRates_SS_Moriond18.root"),
         2018: os.path.join(CJLSTscriptsrelpath, "FakeRates_SS_Moriond19.root"),
      }[year]
    )

__fakerates = KeyDefaultDict(FakeRates)

def getfakerate(year, leppt, lepeta, leplepid):
    return __fakerates[year].GetFakeRate(leppt, lepeta, leplepid)

from ROOT import CRZLLss, test_bit



class normalizeZX(object):
  el = 11**2
  mu = 13**2

  cb_SS = {
    (2016, -mu*mu): 0.9555,
    (2016, -el*el): 1.082,
    (2016, -el*mu): 1.0792,
    (2017, -mu*mu): 1.009,
    (2017, -el*el): 1.379,
    (2017, -el*mu): 1.1087,
    (2018, -mu*mu): 0.9807,
    (2018, -el*el): 1.206,
    (2018, -el*mu): 1.0677,
  }

  fs_ROS_SS = {
    (-el, el): 1.00868,
    (-mu, mu): 1.04015,
    (-el, mu): 1.00823,
    (-mu, el): 1.0049,
  }

  def __new__(cls, year, Z1Flav, Z2Flav):
    return cls.cb_SS[year, Z1Flav*Z2Flav] * cls.fs_ROS_SS[Z1Flav, Z2Flav]
