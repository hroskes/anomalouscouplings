import os

import ROOT

from utilities import KeyDefaultDict
from CJLSTscripts import CJLSTscriptsfolder, convertTGraphstoTH1Fs

class FakeRates(ROOT.FakeRates):
  def __init__(self, arg):
    return super(FakeRates, self).__init__(
      {
         (2016, False): os.path.join(CJLSTscriptsfolder, "FakeRate_SS_Moriond368.root"),
         (2017, False): os.path.join(CJLSTscriptsfolder, "FakeRates_SS_Moriond18.root"),
         (2018, False): os.path.join(CJLSTscriptsfolder, "FakeRates_SS_Moriond19.root"),
         (2016, True): os.path.join(CJLSTscriptsfolder, "FakeRates_SS_2016_Legacy.root"),
         (2017, True): os.path.join(CJLSTscriptsfolder, "FakeRates_SS_2017_Legacy.root"),
         (2018, True): os.path.join(CJLSTscriptsfolder, "FakeRates_SS_2018_Legacy.root"),
      }[arg]
    )

__fakerates = KeyDefaultDict(FakeRates)

def getfakerate(year, uselegacyobjects, leppt, lepeta, leplepid):
    return __fakerates[year, uselegacyobjects].GetFakeRate(leppt, lepeta, leplepid)

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

  ratio_combination_over_SS_newobjects = {
    (2016, -mu*mu): 15.1192528623654 / 13.5568,
    (2016, -el*el): 15.8741490531714 / 13.0417,
    (2016, -el*mu): 33.6281681155334 / 27.8966,
    (2017, -mu*mu): 20.1216409539999 / 17.4629,
    (2017, -el*el): 12.7202599310126 / 10.8171,
    (2017, -el*mu): 32.2092377174772 / 28.2197,
    (2018, -mu*mu): 36.0645009910435 / 33.7859,
    (2018, -el*el): 19.5751513711866 / 16.215,
    (2018, -el*mu): 53.8482753544511 / 47.9655,
  }

  def __new__(cls, year, uselegacyobjects, Z1Flav, Z2Flav):
    if uselegacyobjects: return cls.ratio_combination_over_SS_newobjects[year, Z1Flav*Z2Flav]
    return cls.cb_SS[year, Z1Flav*Z2Flav] * cls.fs_ROS_SS[Z1Flav, Z2Flav]
