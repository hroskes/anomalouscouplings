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
  analysis = Analysis(analysis)
  production = Production(production)
  productionmode = ProductionMode(productionmode)
  if analysis == "fL1Zg" and "L1Zg" in disc and "HadWH" in disc: return

  fainame = analysis.title(superscript=productionmode if productionmode != "ggH" else "dec")
  hff = "Hff0+" if productionmode == "ttH" else None

  hs = {}
  hstack = ROOT.THStack(disc, "")
  l = ROOT.TLegend(.6, .6, .9, .9)

  l.SetBorderSize(0)
  l.SetFillStyle(0)
  stylefunctions.applylegendstyle(l)
  c = pc.TCanvas("c", "", 8, 30, 800, 800)
  stylefunctions.applycanvasstyle(c)

  xaxislabel = discriminant(disc).title

  for color, systematic in enumerate(("Nominal", "ScaleUp", "ScaleDn", "TuneUp", "TuneDn"), start=1):
    if production.year == 2017 and "Scale" in systematic: continue
    sample = ReweightingSamplePlus(productionmode, "0+", hff, "POWHEG", None if systematic == "Nominal" else systematic)
    sample, _, _, reweightfrom = Line(sample, None, 1, production=production)
    if color >= 3: color += 1
    if color >= 5: color += 1
    if color == 7: color = ROOT.kGreen+3
    h = hs[analysis,disc,systematic,sample] = plotfromtree(
      reweightfrom=reweightfrom,
      reweightto=sample,
      disc=disc,
      xaxislabel=xaxislabel,
      normalizeto1=True,
      color=color,
      category=category,
      categorization=analysis.categoryname,
      production=production,
    )
    h.SetLineWidth(2)
    hstack.Add(h)
    l.AddEntry(h, systematic, "l")

  hstack.Draw("hist nostack")
  stylefunctions.applyaxesstyle(hstack)
  l.Draw()
  c.Update()

  hstack.GetXaxis().SetTitle(h.GetXaxis().GetTitle())

  saveasdir = os.path.join(config.plotsbasedir, "xchecks", "pythiasystematics", str(production.year), str(analysis), str(productionmode))
  try:
    os.makedirs(saveasdir)
  except OSError:
    pass
  for ext in "png eps root pdf".split():
    c.SaveAs(os.path.join(saveasdir, "{}.{}".format(disc, ext)))

if __name__ == "__main__":
  with PlotCopier() as pc:
    for analysis in analyses:
      for category in "VBFtagged", "VHHadrtagged":
        for disc in TemplatesFile("ggh", "2e2mu", analysis, category, "180722").discriminants:
          for p in "ggH", "VBF", "ZH", "WplusH", "WminusH", "ttH":
            if p == "ZH": continue
            for production in config.productionsforcombine:
              makeplot(p, analysis, category, production, disc.identifier)
