#!/usr/bin/env python

import os

from helperstuff import config
from helperstuff.copyplots import copyplots
from helperstuff.plotlimits import plotlimits
from helperstuff.utilities import cd, mkdir_p, TFile

def fa3():
  mkdir_p(os.path.join(config.plotsbasedir, "limits", "forfreezing"))
  with cd(os.path.join(config.repositorydir, "scans/fromUlascan/fa3_combined")):
    plotlimits(
      os.path.join(config.plotsbasedir, "limits", "forfreezing", "fa3.png"),
      "fa3",
      0,
      productions=None,
      luminosity=77.3,
      scanranges=((101,-1,1), (100,-0.02,0.02)),
      CMStext="Work in progress",
    )

def fa2():
  with cd(os.path.join(config.repositorydir, "scans/fromUlascan/fa2")):
    plotlimits(
      os.path.join(config.plotsbasedir, "limits", "forfreezing", "fa2.png"),
      "fa2",
      0,
      productions=None,
      luminosity=77.3,
      scanranges=((101,-1,1), (100,-0.02,0.02)),
      CMStext="Work in progress",
      adddirectories=[("/work-zfs/lhc/heshy/bkpanomalouscouplings/HIG-17-011/CMSSW_8_1_0/src/HiggsAnalysis/HZZ4l_Combination/CreateDatacards/FINAL_FOR_HIG-17-011/cards_fa2_allsysts", "_lumi35.8671_13", 35.8671/(2.7+35.8671))],
    )

def fL1():
  with cd(os.path.join("/work-zfs/lhc/heshy/bkpanomalouscouplings/HIG-17-011/CMSSW_8_1_0/src/HiggsAnalysis/HZZ4l_Combination/CreateDatacards/FINAL_FOR_HIG-17-011/cards_fL1_allsysts")):
    plotlimits(
      os.path.join(config.plotsbasedir, "limits", "forfreezing", "fL1.png"),
      "fL1",
      0,
      productions=None,
      luminosity=77.3,
      scanranges=((101,-1,1), (100,-0.02,0.02)),
      CMStext="Work in progress",
      scale=77.3/(2.7+35.8671),
      moreappend="_lumi35.8671_13",
    )

def fL1Zg():
  with cd(os.path.join("/work-zfs/lhc/heshy/bkpanomalouscouplings/HIG-17-011/CMSSW_8_1_0/src/HiggsAnalysis/HZZ4l_Combination/CreateDatacards/FINAL_FOR_HIG-17-011/cards_fL1Zg_allsysts")):
    plotlimits(
      os.path.join(config.plotsbasedir, "limits", "forfreezing", "fL1Zg.png"),
      "fL1Zg",
      0,
      productions=None,
      luminosity=77.3,
      scanranges=((101,-1,1), (100,-0.02,0.02)),
      CMStext="Work in progress",
      scale=77.3/(2.7+35.8671),
      moreappend="_lumi35.8671_13",
    )

if __name__ == "__main__":
  fa3()
  fa2()
  fL1()
  fL1Zg()
  copyplots("limits/forfreezing")
