import os

import ROOT

from utilities import cache, KeyDefaultDict
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

  @classmethod
  @cache
  def ratio_combination_over_SS_newobjects(cls, year, flavor):
    lines = iter(cls.txtfile.split("\n"))
    for line in lines:
      if line.startswith(str(year)): break
    for line in lines:
      if line.startswith("201"): assert False
      if line.startswith({-el*el: "4e", -mu*mu: "4mu", -el*mu: "2e2mu_Combined"}[flavor]):
        SS = float(line.split()[1])
        break
    for line in lines:
      if line.startswith("201"): assert False
      if line.startswith("SS-OS Combination"): break
    for line in lines:
      if line.startswith("201"): assert False
      if line.startswith({-el*el: "4e", -mu*mu: "4mu", -el*mu: "2e2mu"}[flavor]):
        combination = float(line.split()[1])
        break
    return combination / SS

  def __new__(cls, year, uselegacyobjects, Z1Flav, Z2Flav):
    if uselegacyobjects: return cls.ratio_combination_over_SS_newobjects(year, Z1Flav*Z2Flav)
    return cls.cb_SS[year, Z1Flav*Z2Flav] * cls.fs_ROS_SS[Z1Flav, Z2Flav]

  txtfile = """
******************************************************************************
2016		SS		OS	
4mu		26.5 +/- 8.1	25.7 +/- 8.2
4e		13.1 +/- 5.5	20.2 +/- 6.2
2e2mu		21.6 +/- 6.6	23.6 +/- 7.5
2mu2e		16.8 +/- 7.0	23.5 +/- 7.2
2e2mu_Combined	38.5 +/- 9.7	47.1 +/- 10.3
				
SS-OS Combination	Yield	
4mu			26.1	
4e			16.2	
2e2mu			42.5	
				
******************************************************************************
2017		SS		OS	
4mu		30.0 +/- 9.2	31.9 +/- 10.0
4e		10.8 +/- 4.1	16.1 +/- 4.9
2e2mu		23.2 +/- 7.1	22.3 +/- 7.1
2mu2e		14.7 +/- 5.5	20.9 +/- 6.4
2e2mu_Combined 	37.9 +/- 9.0	43.1 +/- 9.5
				
SS-OS Combination	Yield	
4mu			30.9	
4e			13.0	
2e2mu			40.4	

******************************************************************************
2018		SS		OS		
4mu		49.9 +/- 15.2	47.7 +/- 14.8
4e		16.2 +/- 5.9	25.4 +/- 7.7
2e2mu		35.9 +/- 10.9	33.7 +/- 10.6
2mu2e		23.6 +/- 8.5	33.5 +/- 10.1
2e2mu_Combined	59.5 +/- 13.8	67.2 +/- 14.6
				
SS-OS Combination	Yield	
4mu			48.8	
4e			19.6	
2e2mu			63.1	
******************************************************************************
  """
