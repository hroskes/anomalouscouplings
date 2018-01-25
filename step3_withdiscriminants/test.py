#!/usr/bin/env python

assert __name__ == "__main__"

import argparse

p = argparse.ArgumentParser()
p.add_argument("process", choices="decay VBFdecay HadVHdecay".split())
args = p.parse_args()
process = args.process

import glob
import itertools
import os

import ROOT

from helperstuff import config, style
from helperstuff.copyplots import copyplots
from helperstuff.discriminants import discriminant
from helperstuff.treewrapper import TreeWrapper

disc = "D_4couplings_"+process
disc = discriminant(disc)

t = ROOT.TChain("candTree")
for filename in glob.iglob("*.root"):
  if not os.path.exists(filename+".tmp"):
    t.Add(filename)

c = ROOT.TCanvas()
c.SetLogy()
t.Draw(disc.name+">>h({},{},{})".format(disc.bins, disc.min, disc.max), disc.name.replace("D_4couplings", "D_0minus")+">-998")
h = ROOT.h

numbersofbins = [(variablename, len(binning)+1) for variablename, binning in getattr(TreeWrapper, "binning_4couplings_"+process)]
for binnumber in range(disc.bins):
  for D_0minus, D_0hplus, D_L1, D_L1Zg, D_int in itertools.product(
      (0, .5, 1), (0, .6, 1), (0, .6, 1), (0, .5, 1), (-1, 1)
  ):
    result = 0
    for variablename, binning in getattr(TreeWrapper, "binning_4couplings_"+process):
      result *= len(binning)+1
      variable = globals()[variablename.rsplit("_", 1)[0]]
      for bin in binning:
        if variable > bin:
          result += 1
    if result == binnumber:
      if h.GetBinContent(binnumber+1) < 100:
        print "{:3d} {:4.1f} {:4.1f} {:4.1f} {:4.1f} {:4.1f} {:5d}".format(binnumber, D_0minus, D_0hplus, D_L1, D_L1Zg, D_int, int(h.GetBinContent(binnumber+1)))
      break
  else:
    raise ValueError("No bin #{}??".format(binnumber))

for ext in "png eps root pdf".split():
  c.SaveAs(os.path.join(config.plotsbasedir, "TEST", "D_4couplings_"+process+"."+ext))
copyplots("TEST")
