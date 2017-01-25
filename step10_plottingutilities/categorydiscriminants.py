#!/usr/bin/env python

import os

import ROOT

from helperstuff import config, style
from helperstuff.CJLSTscripts import getDZHhWP, getDWHhWP, getDVBF2jetsWP
from helperstuff.discriminants import discriminant
from helperstuff.enums import analyses, Analysis
from helperstuff.samples import ArbitraryCouplingsSample, ReweightingSample
from helperstuff.plotfromtree import Line, plotfromtree

def makeplot(analysis, disc, additionalcconstant=1, shiftWP=None):
  analysis == Analysis(analysis)
  if "HadZH" in disc:
    productionmode = "ZH"
    WP = getDZHhWP(125, config.useQGTagging)
  elif "HadWH" in disc:
    if analysis == "fL1Zg": return
    productionmode = "WH"
    WP = getDWHhWP(125, config.useQGTagging)
  elif "2jet" in disc:
    productionmode = "VBF"
    WP = getDVBF2jetsWP(125, config.useQGTagging)

  if shiftWP is not None:
    if additionalcconstant != 1:
      raise TypeError("Can't provide both arguments additionalcconstant and shiftWP!")
    if not 0 < shiftWP < 1:
      raise ValueError("Can't shift WP to {}!  Has to be between 0 and 1!".format(shiftWP))
    additionalcconstant = WP / shiftWP * (shiftWP-1) / (WP-1)

  WP = (additionalcconstant*(1/WP-1)+1)**-1
  if shiftWP is not None:
    assert abs(WP - shiftWP) < 1e-10, "{} {}".format(WP, shiftWP)

  fainame = "{}^{{{}}}".format(analysis.title, productionmode)

  SM, BSM = (ReweightingSample(productionmode, _) for _ in analysis.purehypotheses)
  mixplus = ReweightingSample(productionmode, analysis.mixprodhypothesis)
  mixminus = ArbitraryCouplingsSample(productionmode, mixplus.g1, -mixplus.g2, -mixplus.g4, -mixplus.g1prime2, -mixplus.ghzgs1prime2)

  lines = [
           Line(SM, "{} SM".format(productionmode), 1, bkpreweightfrom=ReweightingSample(productionmode, "0+")),
           Line(BSM,                            "{} {}=1".format(productionmode, fainame), ROOT.kViolet+7, bkpreweightfrom=ReweightingSample(productionmode, "L1")),
           Line(mixplus,                        "{} {}=#plus0.5".format(productionmode, fainame), ROOT.kGreen+3, bkpreweightfrom=ReweightingSample(productionmode, "L1")),
           Line(mixminus,                       "{} {}=#minus0.5".format(productionmode, fainame), 4, mixplus, bkpreweightfrom=ReweightingSample(productionmode, "L1")),
           Line(ReweightingSample("ggH", "0+"), "ggH", 2),
          ]
  if analysis == "fL1Zg": del lines[3]

  hs = {}
  hstack = ROOT.THStack(disc, "")
  l = ROOT.TLegend(.4, .6, .9, .9)
  l.SetBorderSize(0)
  l.SetFillStyle(0)
  c = ROOT.TCanvas()

  xaxislabel = discriminant(disc).title
  if additionalcconstant != 1:
    xaxislabel += " (shifted)"

  for sample, title, color, reweightfrom in lines:
    h = hs[analysis,disc,sample] = plotfromtree(
      reweightfrom=reweightfrom,
      reweightto=sample,
      disc=disc,
      transformation="1/({cprime}*(1/{{disc}}-1)+1)".format(cprime=additionalcconstant),
      xaxislabel=xaxislabel,
      normalizeto1=True,
      color=color,
    )
    hstack.Add(h)
    l.AddEntry(h, title, "l")

  hstack.Draw("hist nostack")
  l.Draw()

  hstack.GetXaxis().SetTitle(h.GetXaxis().GetTitle())

  line = ROOT.TLine(WP, ROOT.gPad.GetUymin(), WP, l.GetY1())
  line.SetLineWidth(4)
  line.SetLineColor(ROOT.kOrange+7)
  line.SetLineStyle(9)
  line.Draw()
  l.AddEntry(line, "cut", "l")

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
        if (
               analysis == "fa3" and h in ("a2", "L1")
            or analysis == "fa2" and h in ("0minus", "L1")
            or analysis == "fL1" and h in ("0minus", "a2")
           ): continue

        shiftWP = None
        if "Had" in p:
          shiftWP = .5

        makeplot(analysis, "D_{}_{}".format(p, h), shiftWP=shiftWP)
