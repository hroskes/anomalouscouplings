#!/usr/bin/env python

from collections import namedtuple
import os

import ROOT

from helperstuff import config, style
from helperstuff.discriminants import discriminant
from helperstuff.enums import analyses, Analysis
from helperstuff.samples import ArbitraryCouplingsSample, ReweightingSample, Sample
from helperstuff.plotfromtree import plotfromtree

class Line(namedtuple("Line", "sample title color reweightfrom")):
    def __new__(cls, sample, title, color, reweightfrom=None):
        if reweightfrom is None: reweightfrom = sample
        if not isinstance(reweightfrom, Sample):
            assert len(config.productionsforcombine) == 1
            reweightfrom = Sample(reweightfrom, config.productionsforcombine[0])
        return super(Line, cls).__new__(cls, sample, title, color, reweightfrom)

def makeplot(analysis, disc):
  analysis == Analysis(analysis)
  if "HadZH" in disc:
    productionmode = "ZH"
  elif "HadWH" in disc:
    productionmode = "WH"
  elif "2jet" in disc:
    productionmode = "VBF"

  fainame = "{}^{{{}}}".format(analysis.title, productionmode)

  mixplus = ReweightingSample(productionmode, analysis.mixprodhypothesis)
  mixminus = ArbitraryCouplingsSample(productionmode, mixplus.g1, -mixplus.g2, -mixplus.g4, -mixplus.g1prime2)
  print mixminus

  lines = [
           Line(ReweightingSample(productionmode, "0+"), "{} SM".format(productionmode), 1),
           Line(ReweightingSample(productionmode, "0-"), "{} {}=1".format(productionmode, fainame), ROOT.kViolet+7),
           Line(mixplus,                                 "{} {}=#plus0.5".format(productionmode, fainame), ROOT.kGreen+3),
           Line(mixminus,                                "{} {}=#minus0.5".format(productionmode, fainame), 4, mixplus),
           Line(ReweightingSample("ggH", "0+"),          "ggH", 2),
          ]

  hs = {}
  hstack = ROOT.THStack(disc, "")
  l = ROOT.TLegend(.4, .6, .9, .9)
  l.SetBorderSize(0)
  l.SetFillStyle(0)
  c = ROOT.TCanvas()

  for sample, title, color, reweightfrom in lines:
    h = hs[analysis,disc,sample] = plotfromtree(
      productionmode=reweightfrom.productionmode,
      hypothesis=reweightfrom.hypothesis,
      weight=sample,
      disc=disc,
      normalizeto1=True,
      color=color,
    )
    hstack.Add(h)
    l.AddEntry(h, title, "l")

  hstack.Draw("hist nostack")
  hstack.GetXaxis().SetTitle(h.GetXaxis().GetTitle())
  l.Draw()
  saveasdir = os.path.join(config.plotsbasedir, "categorization", str(analysis))
  try:
    os.makedirs(saveasdir)
  except OSError:
    pass
  for ext in "png eps root pdf".split():
    c.SaveAs(os.path.join(saveasdir, "{}.{}".format(disc, ext)))

if __name__ == "__main__":
  for analysis in analyses:
    for p in "HadWH", "HadZH", "2jet":
      for h in "0plus", "0minus", "a2", "L1":
        makeplot(analysis, "D_{}_{}".format(p, h))
