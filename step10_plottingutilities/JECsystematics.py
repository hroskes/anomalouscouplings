#!/usr/bin/env python

import os

import ROOT

from helperstuff import config, stylefunctions
from helperstuff.CJLSTscripts import getDZHhWP, getDWHhWP, getDVBF2jetsWP
from helperstuff.discriminants import discriminant
from helperstuff.enums import analyses, Analysis, Production, ProductionMode
from helperstuff.plotfromtree import Line, plotfromtree
from helperstuff.samples import ArbitraryCouplingsSample, ReweightingSample, ReweightingSamplePlus
from helperstuff.templates import TemplatesFile
from helperstuff.utilities import PlotCopier

def makeplot(productionmode, analysis, category, production, disc):
  if productionmode == "VBF" and category == "VHLepttagged": return
  print productionmode, category, production, disc

  analysis = Analysis(analysis)
  production = Production(production)
  productionmode = ProductionMode(productionmode)

  fainame = analysis.title(superscript=productionmode if productionmode != "ggH" else "dec")
  hff = "Hff0+" if productionmode == "ttH" else None

  if productionmode.issignal:
    hypothesis = "0+"
    SM, _, _, _, _ = (ReweightingSample(productionmode, _, hff) for _ in analysis.purehypotheses)
    if hypothesis: title = str(productionmode) + " SM"
  else:
    title = str(productionmode)
    hypothesis = None
    SM = ReweightingSample(productionmode)

  bkpreweightfrom = None
  if productionmode in ("VBF", "ZH", "WH"): bkpreweightfrom = ReweightingSample(productionmode, "fa3prod0.5", hff)
  if productionmode == "qqZZ": bkpreweightfrom = ReweightingSamplePlus(productionmode, "ext1")
  sample, _, _, reweightfrom = Line(SM, title, 1, bkpreweightfrom=bkpreweightfrom, production=production)

  hs = {}
  hstack = ROOT.THStack(disc, "")
  l = ROOT.TLegend(.6, .6, .9, .9)

  l.SetBorderSize(0)
  l.SetFillStyle(0)
  stylefunctions.applylegendstyle(l)
  c = pc.TCanvas("c", "", 8, 30, 800, 800)
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
      categorization=analysis.categoryname+JECappend,
      production=production,
    )
    h.SetLineWidth(2)
    hstack.Add(h)
    l.AddEntry(h, JEC, "l")

  hstack.Draw("hist nostack")
  stylefunctions.applyaxesstyle(hstack)
  l.Draw()
  c.Update()

  hstack.GetXaxis().SetTitle(h.GetXaxis().GetTitle())

  saveasdir = os.path.join(config.plotsbasedir, "xchecks", "JECsystematics", str(production.year), str(analysis), str(productionmode))
  try:
    os.makedirs(saveasdir)
  except OSError:
    pass
  for ext in "png eps root pdf".split():
    c.SaveAs(os.path.join(saveasdir, "{}.{}".format(disc, ext)))

if __name__ == "__main__":
  with PlotCopier() as pc:
    for analysis in "fa3fa2fL1fL1Zg_morecategories",:
      for category in "VBFtagged", "VHHadrtagged":
        for disc in TemplatesFile("ggh", "2e2mu", analysis, category, "190821_2016").discriminants:
          for p in "ggH", "VBF", "ZH", "WH", "ttH", "qqZZ":
            for production in config.productionsforcombine:
              makeplot(p, analysis, category, production, disc.identifier)
