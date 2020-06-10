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
    lines = iter(cls.txtfilenew.split("\n"))
    for line in lines:
      if line.startswith("# "+str(year)): break
    for line in lines:
      if line.startswith("201"): assert False
      if line.startswith("M4l = [118, 130] GeV"): break
    for line in lines:
      if line.startswith("201"): assert False
      if line.startswith("M4l"): assert False
      if line.startswith({-cls.el*cls.el: "4e", -cls.mu*cls.mu: "4mu", -cls.el*cls.mu: "2e2emu_MERGED"}[flavor]):
        SS = float(line.split()[1])
        combined = float(line.split()[3])
        break
    return combined / SS

  @classmethod
  @cache
  def OS_SS_ratio_new(cls, year, Z1Flav, Z2Flav):
    lines = iter(cls.txtfilenew.split("\n"))
    for line in lines:
      if line.startswith("# "+str(year)): break
    for line in lines:
      if line.startswith("201"): assert False
      if line.startswith("SS/OS ratios"): break  #Matteo: this is a typo, the numbers really are OS/SS
    for line in lines:
      if line.startswith("201"): assert False
      if line.startswith({
        (-cls.el, cls.el): "4e",
        (-cls.el, cls.mu): "2e2mu",
        (-cls.mu, cls.el): "2mu2e",
        (-cls.mu, cls.mu): "4mu",
      }[Z1Flav, Z2Flav]):
        return float(line.split()[1])
    assert False

  def __new__(cls, year, usenewobjects, Z1Flav, Z2Flav):
    if usenewobjects:
      return cls.ratio_combination_over_SS_new(year, Z1Flav*Z2Flav) * cls.OS_SS_ratio_new(year, Z1Flav, Z2Flav)
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
########
# 2016 #
########

M4l > 70 GeV		SS	OS	COMBINATION
4mu			29.66 	26.87	28.19
4e			13.03	20.06	16.13	
2e2emu_MERGED 		41.45	47.45	44.39
2e2mu			24.79	25.18
2mu2e			16.66	22.28


M4l = [118, 130] GeV	SS	OS	COMBINATION
4mu		 	3.66 	3.22	3.48
4e			1.04	1.70	1.24
2e2emu_MERGED 		4.38	4.04	4.24
2e2mu			3.07	2.52
2mu2e			1.31	1.52


SS/OS ratios
4mu	1.00008 +/- 0.027524
4e    	1.00229 +/- 0.011872
2e2mu 	1.03601 +/- 0.030332
2mu2e 	0.99873 +/- 0.010537

########
# 2017 #
########

M4l > 70 GeV		SS	OS	COMBINATION
4mu			33.55 	32.72	33.13
4e			10.91	16.15	12.95	
2e2emu_MERGED 		40.96	45.29	43.05
2e2mu			26.27	23.97
2mu2e			14.69	21.32


M4l = [118, 130] GeV	SS	OS	COMBINATION
4mu		 	4.23 	3.54	3.92
4e			1.02	1.54	1.18
2e2emu_MERGED 		4.55	4.92	4.69
2e2mu			3.20	3.03
2mu2e			1.35	1.90


SS/OS ratios
4mu   1.03961 +/- 0.0267732
4e    1.01168 +/- 0.0133765
2e2mu 1.01307 +/- 0.0301085
2mu2e 1.00266 +/- 0.0115351

########
# 2018 #
########

M4l > 70 GeV		SS	OS	COMBINATION
4mu	 		52.17 	49.38	50.72
4e			15.99	25.35	19.42	
2e2emu_MERGED 		60.78	67.14	63.87
2e2mu			37.53	34.18
2mu2e			23.25	32.96


M4l = [118, 130] GeV	SS	OS	COMBINATION
4mu		 	6.54 	5.59	6.08
4e			1.43	2.16	1.65
2e2emu_MERGED 		6.93	6.77	6.86
2e2mu			4.83	3.99
2mu2e			2.10	2.78


SS/OS ratios
4mu	1.02763 +/- 0.021062
4e    	1.00635 +/- 0.011001
2e2mu 	1.03170 +/- 0.024905
2mu2e 	1.00492 +/- 0.009098
  """
