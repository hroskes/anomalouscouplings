#!/usr/bin/env python

import os

import ROOT

from helperstuff import config, style
from helperstuff.CJLSTscripts import getDZHhWP, getDWHhWP, getDVBF2jetsWP
from helperstuff.discriminants import discriminant
from helperstuff.enums import analyses, Analysis
from helperstuff.samples import ArbitraryCouplingsSample, ReweightingSample
from helperstuff.plotfromtree import Line, plotfromtree

def makeplot(analysis, disc):
  analysis == Analysis(analysis)
  WP = 0.5
  if "HadZH" in disc:
    productionmode = "ZH"
  elif "HadWH" in disc:
    if analysis == "fL1Zg" and "L1Zg" in disc: return
    productionmode = "WH"
  elif "2jet" in disc:
    productionmode = "VBF"

  fainame = "{}^{{{}}}".format(analysis.title, productionmode)

  SM, BSM = (ReweightingSample(productionmode, _) for _ in analysis.purehypotheses)
  mixplus = ReweightingSample(productionmode, analysis.mixprodhypothesis)
  mixminus = ReweightingSample(productionmode, str(analysis.mixprodhypothesis).replace("0.5", "-0.5"))

  lines = [
           Line(SM, "{} SM".format(productionmode), 1, bkpreweightfrom=ReweightingSample(productionmode, "0+")),
           Line(BSM,                            "{} {}=1".format(productionmode, fainame), ROOT.kViolet+7, bkpreweightfrom=ReweightingSample(productionmode, "L1")),
           Line(mixplus,                        "{} {}=#plus0.5".format(productionmode, fainame), ROOT.kGreen+3, bkpreweightfrom=ReweightingSample(productionmode, "L1")),
           Line(mixminus,                       "{} {}=#minus0.5".format(productionmode, fainame), 4, mixplus, bkpreweightfrom=ReweightingSample(productionmode, "L1")),
           Line(ReweightingSample("ggH", "0+"), "ggH", 2),
          ]

  if productionmode == "WH" and analysis == "fL1Zg": del lines[1:4]

  hs = {}
  hstack = ROOT.THStack(disc, "")
  l = ROOT.TLegend(.4, .6, .9, .9)
  l.SetBorderSize(0)
  l.SetFillStyle(0)
  c = ROOT.TCanvas()

  xaxislabel = discriminant(disc).title

  for sample, title, color, reweightfrom in lines:
    h = hs[analysis,disc,sample] = plotfromtree(
      reweightfrom=reweightfrom,
      reweightto=sample,
      disc=disc,
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
      for h in "0plus", "0minus", "a2", "L1", "L1Zg":
        if h.replace("plus", "+").replace("minus", "-") not in analysis.purehypotheses and h != "0plus":
           continue

        makeplot(analysis, "D_{}_{}".format(p, h))
