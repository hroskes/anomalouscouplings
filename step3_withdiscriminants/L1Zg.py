#!/usr/bin/env python

from collections import namedtuple
import os

import ROOT

from helperstuff import config, style
from helperstuff.CJLSTscripts import getDZHhWP, getDWHhWP, getDVBF2jetsWP
from helperstuff.discriminants import discriminant
from helperstuff.enums import analyses, Analysis, ProductionMode
from helperstuff.samples import ArbitraryCouplingsSample, ReweightingSample, Sample
from helperstuff.plotfromtree import plotfromtree

class Line(namedtuple("Line", "sample title color reweightfrom")):
    def __new__(cls, sample, title, color, reweightfrom=None):
        if reweightfrom is None: reweightfrom = sample
        if not isinstance(reweightfrom, Sample):
            assert len(config.productionsforcombine) == 1
            reweightfrom = Sample(reweightfrom, config.productionsforcombine[0])
        return super(Line, cls).__new__(cls, sample, title, color, reweightfrom)

def makeplot(productionmode, disc):
  analysis = Analysis("fL1Zg")
  productionmode = ProductionMode(productionmode)
  if productionmode == "ZH":
    fainame = "f_{#Lambda1Z#gamma}^{ZH}"
    mixplus = ReweightingSample(productionmode, analysis.mixprodhypothesis)
  elif productionmode == "VBF":
    fainame = "f_{#Lambda1Z#gamma}^{VBF}"
    mixplus = ReweightingSample(productionmode, analysis.mixprodhypothesis)
  elif productionmode == "ggH":
    fainame = "f_{#Lambda1Z#gamma}^{dec}"
    mixplus = ReweightingSample(productionmode, analysis.mixdecayhypothesis)

  mixminus = ArbitraryCouplingsSample(productionmode, mixplus.g1, -mixplus.g2, -mixplus.g4, -mixplus.g1prime2, -mixplus.ghzgs1prime2)

  SM = ReweightingSample(productionmode, "0+")

  lines = [
           Line(ReweightingSample(productionmode, "0+"),   "{} SM".format(productionmode), 1),
           Line(ReweightingSample(productionmode, "L1Zg"), "{} {}=1".format(productionmode, fainame), ROOT.kViolet+7, SM),
           Line(mixplus,                                   "{} {}=#plus0.5".format(productionmode, fainame), ROOT.kGreen+3, SM),
#           Line(mixminus,                                  "{} {}=#minus0.5".format(productionmode, fainame), 4, SM),
          ]
  if productionmode == "ZH":
    lines.append(
                 Line(ReweightingSample("WH", "0+"),            "WH SM", 2),
                )

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
  l.Draw()

  hstack.GetXaxis().SetTitle(h.GetXaxis().GetTitle())

  saveasdir = os.path.join(config.plotsbasedir, "TEST")
  try:
    os.makedirs(saveasdir)
  except OSError:
    pass
  for ext in "png eps root pdf".split():
    c.SaveAs(os.path.join(saveasdir, "{}.{}".format(disc, ext)))

if __name__ == "__main__":
  makeplot("VBF", "D_L1Zg_VBFdecay")
  makeplot("VBF", "D_L1Zg_VBF")
  makeplot("VBF", "D_L1Zgint_VBF")
  makeplot("ZH", "D_L1Zg_HadVH")
  makeplot("ZH", "D_L1Zg_HadVHdecay")
  makeplot("ZH", "D_L1Zgint_HadVH")
  makeplot("ggH", "D_L1Zg_decay")
  makeplot("ggH", "D_L1Zgint_decay")
