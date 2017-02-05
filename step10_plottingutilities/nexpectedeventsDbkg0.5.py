#!/usr/bin/env python

import os

import ROOT

from helperstuff import config
from helperstuff.enums import analyses, channels, categories
from projections import Projections

assert len(config.productionsforcombine) == 1
production = config.productionsforcombine[0]

class NExpected(Projections):
  @property
  def nexpected(self):
    disc = self.discriminants[0]
    rootfile = os.path.join(self.saveasdir, disc.name+".root")
    f = ROOT.TFile(rootfile)
    c = f.c1
    hstack = c.GetListOfPrimitives()[1]
    total = 0
    for h in hstack.GetHists():
      if h.GetLineColor() in (1, 6, 2, ROOT.kOrange+6, ROOT.kViolet+7): total += h.Integral()
      elif h.GetLineColor() in (ROOT.kCyan, ROOT.kGreen+3, 4): pass
      else: assert False
    return total

if __name__ == "__main__":
  for analysis in analyses:
    print sum(NExpected(production, channel, category, "rescalemixtures", analysis, "fullrange").nexpected for channel in channels for category in categories),
  print
  for analysis in analyses:
    print sum(NExpected(production, channel, category, "rescalemixtures", analysis, "enrich").nexpected for channel in channels for category in categories),
