#!/usr/bin/env python

import glob, os, shutil, tempfile

import numpy as np

from helperstuff.utilities import TFile, KeepWhileOpenFile

def run(filename):
  oldfilename = filename
  newfilename = filename.replace("bkp/bkp_190821_badweight/", "")
  print newfilename
  assert newfilename != oldfilename
  with KeepWhileOpenFile(newfilename+".tmp") as kwof, tempfile.NamedTemporaryFile(suffix=".root", dir="/home-2/jroskes1@jhu.edu/scratch/tmparea/") as tempf:
    if not kwof: return
    if os.path.exists(newfilename): return
    if "ZX" in newfilename:
      with TFile(oldfilename) as oldf:
        t = oldf.candTree
        t.GetEntry(0)
        assert not hasattr(t, "L1prefiringWeight")
        shutil.copy(oldfilename, newfilename)
        return
    os.remove(tempf.name)
    try:
      with TFile(oldfilename) as oldf, TFile(tempf.name, "CREATE", deleteifbad=True) as newf:
        t = oldf.candTree
        newt = t.CloneTree(0, "fast")
        newt.SetDirectory(newf)
        branch = np.array([0], dtype=np.float64)
        newt.SetBranchAddress("MC_weight_nominal", branch)

        nentries = t.GetEntries()
        for i, entry in enumerate(t, start=1):
          branch[0] = t.MC_weight_nominal * t.L1prefiringWeight
          newt.Fill()
          if i%10000 == 0:
            print i, "/", nentries
      shutil.move(tempf.name, newfilename)
    finally:
      with open(tempf.name, "w"): pass

for filename in sorted(glob.iglob("bkp/bkp_190821_badweight/*/*.root"), key=lambda x: ("qqZZ" not in x) - 2*("4tau" in x)):
  run(filename)
