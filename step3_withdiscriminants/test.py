#!/usr/bin/env python

import glob
import itertools
import os

import ROOT

from helperstuff import config, style
from helperstuff.copyplots import copyplots
from helperstuff.discriminants import discriminant
from helperstuff.treewrapper import TreeWrapper

disc = "D_4couplings_decay"

assert __name__ == "__main__"

disc = discriminant(disc)

t = ROOT.TChain("candTree")
for filename in glob.iglob("*.root"):
  if not os.path.exists(filename+".tmp"):
    t.Add(filename)

c = ROOT.TCanvas()
c.SetLogy()
t.Draw(disc.name+">>h({},{},{})".format(disc.bins, disc.min, disc.max))
h = ROOT.h

numbersofbins = [(variablename, len(binning)+1) for variablename, binning in TreeWrapper.binning_4couplings_decay]
for binnumber in range(disc.bins):
  for D_0minus_decay, D_0hplus_decay, D_L1_decay, D_L1Zg_decay, D_int_decay in itertools.product(
      (0, .5, 1),     (0, .6, 1),     (0, .6, 1), (0, .5, 1),   (-1, 1)
  ):
    result = 0
    for variablename, binning in TreeWrapper.binning_4couplings_decay:
      result *= len(binning)+1
      variable = globals()[variablename]
      for bin in binning:
        if variable > bin:
          result += 1
    if result == binnumber:
      if h.GetBinContent(binnumber+1) < 100:
        print "{:3d} {:4.1f} {:4.1f} {:4.1f} {:4.1f} {:4.1f} {:5d}".format(binnumber, D_0minus_decay, D_0hplus_decay, D_L1_decay, D_L1Zg_decay, D_int_decay, int(h.GetBinContent(binnumber+1)))
      break
  else:
    raise ValueError("No bin #{}??".format(binnumber))

for ext in "png eps root pdf".split():
  c.SaveAs(os.path.join(config.plotsbasedir, "TEST", "test."+ext))
copyplots("TEST")
