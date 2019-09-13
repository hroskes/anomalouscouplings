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
            print otherfa3, otherfa2, otherfL1, otherfL1Zg

            setparametersforgrid = []
            whatsleft = 1
            if otherfa3:
              otherfa3rel = float("{:.1e}".format(otherfa3 / whatsleft))
              if otherfa3rel: otherfa3 = otherfa3rel * whatsleft
              whatsleft -= abs(otherfa3)
              print "fa3", otherfa3rel,
              setparametersforgrid.append("CMS_zz4l_fai1_relative:{}".format(otherfa3rel))
            if otherfa2:
              otherfa2rel = float("{:.1e}".format(otherfa2 / whatsleft))
              if otherfa2rel: otherfa2 = otherfa2rel * whatsleft
              whatsleft -= abs(otherfa2)
              print "fa2", otherfa2rel,
              setparametersforgrid.append("CMS_zz4l_fai2_relative:{}".format(otherfa2rel))
            if otherfL1:
              otherfL1rel = float("{:.1e}".format(otherfL1 / whatsleft))
              if otherfL1rel: otherfL1 = otherfL1rel * whatsleft
              whatsleft -= abs(otherfL1)
              print "fL1", otherfL1rel,
              setparametersforgrid.append("CMS_zz4l_fai3_relative:{}".format(otherfL1rel))
            if otherfL1Zg:
              otherfL1Zgrel = float("{:.1e}".format(otherfL1Zg / whatsleft))
              if otherfL1Zgrel: otherfL1Zg = otherfL1Zgrel * whatsleft
              whatsleft -= abs(otherfL1Zg)
              print "fL1Zg", otherfL1Zgrel,
              setparametersforgrid.append("CMS_zz4l_fai4_relative:{}".format(otherfL1Zgrel))

            print

            priority = 0
            if fai != "fa3": priority += otherfa3**2
            if fai != "fa2": priority += min((otherfa2+0.005)**2, (otherfa2 + 0.3)**2)
            if fai != "fL1": priority += min((otherfL1-0.02)**2, (otherfL1 - 0.15)**2)
            if fai != "fL1Zg": priority += min((otherfL1Zg-0.02)**2, (otherfL1Zg + 0.05)**2)
            #print priority

            jobs.append((priority, [
              "sbatch", "--mem", "12G", "--time", "1:0:0", "slurm.sh",
              "./step9_runcombine.py", "fa3fa2fL1fL1Zg_morecategories", "applySIPcut", "scanranges="+scanrange, "plotnuisances="+plotnuisances,
              "scanfai="+fai,
              "setparametersforgrid="+",".join(setparametersforgrid),
              "categories=Untagged", "usesystematics=0", "expectvalues="
            ]))
            #print setparametersforgrid, priority

jobs.sort()
with open("data/alreadyrun", "a") as f: pass
with open("data/alreadyrun") as f:
  alreadyrun = set(f.read().split("\n"))

try:
  for i, (priority, job) in enumerate(jobs):
    if " ".join(job) in alreadyrun: continue
    alreadyrun.add(" ".join(job))
    print job
    subprocess.check_call(job)
finally:
  with open("data/alreadyrun", "a") as f:
    f.write("\n".join(alreadyrun) + "\n")
