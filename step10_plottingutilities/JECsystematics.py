#!/usr/bin/env python

import os

import ROOT

from helperstuff import config, stylefunctions
from helperstuff.CJLSTscripts import getDZHhWP, getDWHhWP, getDVBF2jetsWP
from helperstuff.discriminants import discriminant
from helperstuff.enums import analyses, Analysis
from helperstuff.plotfromtree import Line, plotfromtree
from helperstuff.samples import ArbitraryCouplingsSample, ReweightingSample
from helperstuff.templates import TemplatesFile

def makeplot(productionmode, analysis, category, disc):
  analysis = Analysis(analysis)
  if analysis == "fL1Zg" and "L1Zg" in disc and "HadWH" in disc: return

  fainame = analysis.title(superscript=productionmode if productionmode != "ggH" else "dec")
  hff = "Hff0+" if productionmode == "ttH" else None

  SM, BSM = (ReweightingSample(productionmode, _, hff) for _ in analysis.purehypotheses)

  sample, _, _, reweightfrom = Line(SM, "{} SM".format(productionmode), 1, bkpreweightfrom=ReweightingSample(productionmode, "0+", hff))

  hs = {}
  hstack = ROOT.THStack(disc, "")
  l = ROOT.TLegend(.6, .6, .9, .9)

  l.SetBorderSize(0)
  l.SetFillStyle(0)
  stylefunctions.applylegendstyle(l)
  c = ROOT.TCanvas("c", "", 8, 30, 800, 800)
  stylefunctions.applycanvasstyle(c)

  xaxislabel = discriminant(disc).title

  for color, JEC in enumerate(("Nominal", "JECUp", "JECDn"), start=1):
    if color == 3: color = 4
    JECappend = ("" if JEC == "Nominal" else "_"+JEC)
    h = hs[analysis,disc,JEC,sample] = plotfromtree(
      reweightfrom=reweightfrom,
      reweightto=sample,
      disc=disc+JECappend,
      xaxislabel=xaxislabel,
      normalizeto1=True,
      color=color,
      category=category,
      analysis=analysis
    )
    h.SetLineWidth(2)
    hstack.Add(h)
    l.AddEntry(h, JEC, "l")

  hstack.Draw("hist nostack")
  stylefunctions.applyaxesstyle(hstack)
  l.Draw()
  c.Update()

  hstack.GetXaxis().SetTitle(h.GetXaxis().GetTitle())

  saveasdir = os.path.join(config.plotsbasedir, "xchecks", "JECsystematics", str(analysis), str(productionmode))
  try:
    os.makedirs(saveasdir)
  except OSError:
    pass
  for ext in "png eps root pdf".split():
    c.SaveAs(os.path.join(saveasdir, "{}.{}".format(disc, ext)))

if __name__ == "__main__":
  for analysis in analyses:
    for category in "VBFtagged", "VHHadrtagged":
      for disc in TemplatesFile("ggh", "2e2mu", analysis, category).discriminants[0:2]:
        for p in "ggH", "VBF", "ZH", "WH", "ttH":
          makeplot(p, analysis, category, disc.name)
