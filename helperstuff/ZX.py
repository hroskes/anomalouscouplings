import os

import ROOT

from CJLSTscripts import CJLSTscriptsfolder, convertTGraphstoTH1Fs
from utilities import cache, KeyDefaultDict

class FakeRates(ROOT.FakeRates):
  def __init__(self, arg):
    return super(FakeRates, self).__init__(
      {
         (2016, False): os.path.join(CJLSTscriptsfolder, "FakeRates_SS_2016_Legacy.root"),
         (2017, False): os.path.join(CJLSTscriptsfolder, "FakeRates_SS_2017_Legacy.root"),
         (2018, False): os.path.join(CJLSTscriptsfolder, "FakeRates_SS_2018_Legacy.root"),
         (2016, True): os.path.join(CJLSTscriptsfolder, "newData_FakeRates_SS_2016.root"),
         (2017, True): os.path.join(CJLSTscriptsfolder, "newData_FakeRates_SS_2017.root"),
         (2018, True): os.path.join(CJLSTscriptsfolder, "newData_FakeRates_SS_2018.root"),
      }[arg]
    )

__fakerates = KeyDefaultDict(FakeRates)

def getfakerate(year, usenewobjects, leppt, lepeta, leplepid):
    return __fakerates[year, usenewobjects].GetFakeRate(leppt, lepeta, leplepid)

from ROOT import CRZLLss, test_bit



class normalizeZX(object):
  el = 11**2
  mu = 13**2

  @classmethod
  @cache
  def ratio_combination_over_SS_old(cls, year, flavor):
    lines = iter(cls.txtfileold.split("\n"))
    for line in lines:
      if line.startswith(str(year)): break
    for line in lines:
      if line.startswith("201"): assert False
      if line.startswith({-cls.el*cls.el: "4e", -cls.mu*cls.mu: "4mu", -cls.el*cls.mu: "2e2mu_Combined"}[flavor]):
        SS = float(line.split()[1])
        break
    for line in lines:
      if line.startswith("201"): assert False
      if line.startswith("SS-OS Combination"): break
    for line in lines:
      if line.startswith("201"): assert False
      if line.startswith({-cls.el*cls.el: "4e", -cls.mu*cls.mu: "4mu", -cls.el*cls.mu: "2e2mu"}[flavor]):
        combination = float(line.split()[1])
        break
    return combination / SS

  @classmethod
  @cache
  def ratio_combination_over_SS_new(cls, year, flavor):
    return {
      2016: {
        -cls.mu*cls.mu: 0.9504,
        -cls.el*cls.el: 1.2379,
        -cls.el*cls.mu: 1.0709,
      },
      2017: {
        -cls.mu*cls.mu: 0.9875,
        -cls.el*cls.el: 1.1870,
        -cls.el*cls.mu: 1.0510,
      },
      2018: {
        -cls.mu*cls.mu: 0.9722,
        -cls.el*cls.el: 1.2145,
        -cls.el*cls.mu: 1.05088,
      },
    }[year][flavor]

  def __new__(cls, year, usenewobjects, Z1Flav, Z2Flav):
    if usenewobjects:
      return cls.ratio_combination_over_SS_new(year, Z1Flav*Z2Flav)
    else:
      return cls.ratio_combination_over_SS_old(year, Z1Flav*Z2Flav)

  txtfileold = """
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

  txtfilenew = """
SS	2016	2017	2018	Full Run 2
4mu	9.48	10.94	16.87	37.29		
4e	2.67	2.63	4.10	9.40
2e2mu	8.26	8.18	13.26	29.70
2mu2e	3.46	3.52	5.55	12.53
TOT	23.87	25.27	39.78	88.92

OS	2016	2017	2018	Full Run 2
4mu	8.37	8.58	14.11	31.06		
4e	5.33	3.85	6.53	15.71
2e2mu	7.68	7.59	11.05	26.32
2mu2e	5.67	6.39	8.68	20.74
TOT	27.05	26.41	40.37	93.83

COMB	2016	2017	2018	Full Run 2
4mu	8.94	9.66	15.37	33.97		
4e	3.43	3.07	4.88	11.38
2e2mu	12.43	12.63	19.25	44.31
TOT	24.8	25.36	39.5	89.66
  """
