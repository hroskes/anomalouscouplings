#!/usr/bin/env python

import os, tempfile

from helperstuff import config
from helperstuff.copyplots import copyplots
from helperstuff.plotlimits import plotlimits
from helperstuff.utilities import cd, mkdir_p, TFile

from mergeplots import mergeplots

tmp = tempfile.mkdtemp()

def fa3():
  mkdir_p(os.path.join(config.plotsbasedir, "limits", "forcombJune1"))
  mkdir_p(os.path.join(tmp, "fa3_all"))
  mkdir_p(os.path.join(tmp, "fa3_13TeV"))
  with cd(os.path.join(config.repositorydir, "scans/fromUlascan/fa3_2017")):
    plotlimits(
      os.path.join(tmp, "fa3_all", "fa3.png"),
      "fa3",
      0,
      productions=None,
      luminosity=77.3,
      scanranges=((101,-1,1), (100,-0.02,0.02)),
      CMStext="Work in progress",
      adddirectories=[("../fa3_20152016", "_lumi35.8671_13", 1), "../fa3_Run1"],
    )

    plotlimits(
      os.path.join(tmp, "fa3_13TeV", "fa3.png"),
      "fa3",
      0,
      productions=None,
      luminosity=77.3,
      scanranges=((101,-1,1), (100,-0.02,0.02)),
      CMStext="Work in progress",
      adddirectories=[("../fa3_20152016", "_lumi35.8671_13", 1)],
    )

    mergeplots("fa3", subdir=tmp, lumi=77.3, outdir=os.path.join(config.plotsbasedir, "limits", "forcombJune1"), copy=False)

def fa2():
  mkdir_p(os.path.join(config.plotsbasedir, "limits", "forcombJune1"))
  mkdir_p(os.path.join(tmp, "fa2_all"))
  mkdir_p(os.path.join(tmp, "fa2_13TeV"))
  with cd(os.path.join(config.repositorydir, "scans/fromUlascan/fa2_2017")):
    plotlimits(
      os.path.join(tmp, "fa2_all", "fa2.png"),
      "fa2",
      0,
      productions=None,
      luminosity=77.3,
      scanranges=((101,-1,1), (100,-0.02,0.02)),
      CMStext="Work in progress",
      adddirectories=[("../fa2_20152016", "_lumi35.8671_13", 1), "../fa2_Run1"],
      scale=0,
    )

    plotlimits(
      os.path.join(tmp, "fa2_13TeV", "fa2.png"),
      "fa2",
      0,
      productions=None,
      luminosity=77.3,
      scanranges=((101,-1,1), (100,-0.02,0.02)),
      CMStext="Work in progress",
      adddirectories=[("../fa2_20152016", "_lumi35.8671_13", 1)],
      scale=0,
    )

    mergeplots("fa2", subdir=tmp, lumi=77.3, outdir=os.path.join(config.plotsbasedir, "limits", "forcombJune1"), copy=False)

def fL1():
  mkdir_p(os.path.join(config.plotsbasedir, "limits", "forcombJune1"))
  mkdir_p(os.path.join(tmp, "fL1_all"))
  mkdir_p(os.path.join(tmp, "fL1_13TeV"))
  with cd(os.path.join(config.repositorydir, "scans/fromUlascan/fL1_20152016")):
    plotlimits(
      os.path.join(tmp, "fL1_all", "fL1.png"),
      "fL1",
      0,
      productions=None,
      luminosity=77.3,
      scanranges=((101,-1,1), (100,-0.02,0.02)),
      CMStext="Work in progress",
      scale=77.3/(2.7+35.8671),
      adddirectories=[("../fL1_Run1", "", 1)],
    )

    plotlimits(
      os.path.join(tmp, "fL1_13TeV", "fL1.png"),
      "fL1",
      0,
      productions=None,
      luminosity=77.3,
      scanranges=((101,-1,1), (100,-0.02,0.02)),
      CMStext="Work in progress",
      scale=77.3/(2.7+35.8671),
      adddirectories=[],
    )

    mergeplots("fL1", subdir=tmp, lumi=77.3, outdir=os.path.join(config.plotsbasedir, "limits", "forcombJune1"), copy=False)

def fL1Zg():
  mkdir_p(os.path.join(config.plotsbasedir, "limits", "forcombJune1"))
  mkdir_p(os.path.join(tmp, "fL1Zg_all"))
  mkdir_p(os.path.join(tmp, "fL1Zg_13TeV"))
  with cd(os.path.join(config.repositorydir, "scans/fromUlascan/fL1Zg_20152016")):
    plotlimits(
      os.path.join(tmp, "fL1Zg_all", "fL1Zg.png"),
      "fL1Zg",
      0,
      productions=None,
      luminosity=77.3,
      scanranges=((101,-1,1), (100,-0.02,0.02)),
      CMStext="Work in progress",
      scale=77.3/(2.7+35.8671) - 1,
      adddirectories=[("../fL1Zg_Run120152016", "_lumi35.8671_7813", 1)],
    )

    plotlimits(
      os.path.join(tmp, "fL1Zg_13TeV", "fL1Zg.png"),
      "fL1Zg",
      0,
      productions=None,
      luminosity=77.3,
      scanranges=((101,-1,1), (100,-0.02,0.02)),
      CMStext="Work in progress",
      scale=77.3/(2.7+35.8671),
      adddirectories=[],
    )

    mergeplots("fL1Zg", subdir=tmp, lumi=77.3, outdir=os.path.join(config.plotsbasedir, "limits", "forcombJune1"), copy=False)

if __name__ == "__main__":
  fa3()
  fa2()
  fL1()
  fL1Zg()
  copyplots("limits/forcombJune1")
