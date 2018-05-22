#!/usr/bin/env python

import os

from helperstuff import config
from helperstuff.copyplots import copyplots
from helperstuff.plotlimits import plotlimits
from helperstuff.utilities import cd, mkdir_p

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
  pass

def fL1():
  pass

def fL1Zg():
  pass

if __name__ == "__main__":
  fa3()
  fa2()
  fL1()
  fL1Zg()
  copyplots("limits/forfreezing")
