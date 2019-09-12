#!/usr/bin/env python

import itertools
import numpy as np
import subprocess

otherfa3s   = [[0, 0.1, 0.2]]
otherfa2s   = [np.arange(-1, 1.1, 0.2), np.arange(-.4, -.19, 0.02), np.arange(-0.01, 0, 0.0025)]
otherfL1s   = [np.arange(-1, 1.1, 0.2), np.arange(.1, .21, 0.02), np.arange(0.015, 0.0275, 0.0025)]
otherfL1Zgs = [np.arange(-1, 1.1, 0.2), np.arange(-0.07, -0.02, 0.01), np.arange(0, 0.05, 0.01)]

otherfa3s = np.array(sorted(set(np.round(np.concatenate(otherfa3s), 3))))
otherfa2s = np.array(sorted(set(np.round(np.concatenate(otherfa2s), 3))))
otherfL1s = np.array(sorted(set(np.round(np.concatenate(otherfL1s), 3))))
otherfL1Zgs = np.array(sorted(set(np.round(np.concatenate(otherfL1Zgs), 3))))

otherfais = {"fa3": otherfa3s, "fa2": otherfa2s, "fL1": otherfL1s, "fL1Zg": otherfL1Zgs}
fais = "fa3", "fa2", "fL1", "fL1Zg"
scanranges = "101,-1,1", "101,-0.02,0.02"

jobs = []

for scanrange in scanranges:
  for fai in fais:
    plotnuisances = ["CMS_zz4l_fai1", "CMS_zz4l_fai2", "CMS_zz4l_fai3", "CMS_zz4l_fai4", "CMS_zz4l_fa1"]
    del plotnuisances[fais.index(fai)]
    plotnuisances = ",".join(plotnuisances)
    for otherfa3 in (otherfais["fa3"] if fai != "fa3" else [0]):
      for otherfa2 in (otherfais["fa2"] if fai != "fa2" else [0]):
        for otherfL1 in (otherfais["fL1"] if fai != "fL1" else [0]):
          for otherfL1Zg in (otherfais["fL1Zg"] if fai != "fL1Zg" else [0]):
            if abs(otherfa3) + abs(otherfa2) + abs(otherfL1) + abs(otherfL1Zg) > 1: continue
            #print otherfa3, otherfa2, otherfL1, otherfL1Zg,

            setparametersforgrid = []
            if otherfa3: setparametersforgrid.append("CMS_zz4l_fai1_relative:"+str(otherfa3))
            if otherfa2: setparametersforgrid.append("CMS_zz4l_fai2_relative:"+str(otherfa2))
            if otherfL1: setparametersforgrid.append("CMS_zz4l_fai3_relative:"+str(otherfL1))
            if otherfL1Zg: setparametersforgrid.append("CMS_zz4l_fai4_relative:"+str(otherfL1Zg))

            priority = 0
            if fai != "fa3": priority += otherfa3**2
            if fai != "fa2": priority += min((otherfa2+0.005)**2, (otherfa2 + 0.3)**2)
            if fai != "fL1": priority += min((otherfL1-0.02)**2, (otherfL1 - 0.15)**2)
            if fai != "fL1Zg": priority += min((otherfL1Zg-0.02)**2, (otherfL1Zg + 0.05)**2)
            #print priority

            jobs.append((priority, [
              "sbatch", "--mem", "12G", "--time", "6:0:0", "slurm.sh",
              "./step9_runcombine.py", "fa3fa2fL1fL1Zg_morecategories", "applySIPcut", "scanranges="+scanrange, "plotnuisances="+plotnuisances,
              "useNLLandNLL0=0", "scanfai="+fai,
              "setparametersforgrid="+",".join(setparametersforgrid),
              "categories=Untagged", "usesystematics=0", "expectvalues="
            ]))
            #print setparametersforgrid, priority

jobs.sort()
with open("data/alreadyrun", "a") as f: pass
with open("data/alreadyrun") as f:
  alreadyrun = set(f.read().split("\n"))

try:
  for priority, job in jobs:
    print priority
    if " ".join(job) in alreadyrun: continue
    alreadyrun.add(" ".join(job))
    print job
    subprocess.check_call(job)
finally:
  with open("data/alreadyrun", "a") as f:
    f.write("\n".join(alreadyrun) + "\n")
