#!/usr/bin/env python

import os

import ROOT

from helperstuff import config, style
from helperstuff.enums import analyses, categories
from helperstuff.plotfromtree import Line, plotfromtree
from helperstuff.samples import ReweightingSample
from helperstuff.templates import TemplatesFile
from helperstuff.utilities import mkdir_p

lines = [
         Line(ReweightingSample("qqZZ"), "qq#rightarrowZZ", 1),
         Line(ReweightingSample("VBF bkg", "2e2mu"), "VBF bkg", 2),
         Line(ReweightingSample("VBF bkg", "4e"), "VBF bkg", 2),
         Line(ReweightingSample("qqZZ", "2e2mu"), "VBF bkg rwt to qq#rightarrowZZ", 4, reweightfrom=ReweightingSample("VBF bkg", "2e2mu")),
         Line(ReweightingSample("qqZZ", "4e"), "VBF bkg rwt to qq#rightarrowZZ", 4, reweightfrom=ReweightingSample("VBF bkg", "4e")),
        ]

def plot(disc, category, analysis):
  hs = {}
  hstack = ROOT.THStack(disc, "")
  l = ROOT.TLegend(.4, .5, .85, .9)
  l.SetBorderSize(0)
  l.SetFillStyle(0)
  c = ROOT.TCanvas()

  for sample, title, color, reweightfrom in lines:
    h = hs[analysis,disc,sample] = plotfromtree(
      reweightfrom=reweightfrom,
      reweightto=sample,
      disc=disc,
      normalizeto1=False,
      color=color,
      #channel="2e2mu",
      category=category,
      analysis=analysis,
    )
    if sample.flavor == "4e":
      hs[analysis,disc,ReweightingSample(sample.productionmode, "2e2mu")].Add(h)
      continue
    h.SetLineWidth(2)
    hstack.Add(h)
    l.AddEntry(h, title, "l")

  for h in hs.values():
    h.Scale(1/h.Integral())

  hstack.Draw("hist nostack")
  l.Draw()
  hstack.GetXaxis().SetTitle(h.GetXaxis().GetTitle())

  saveasdir = os.path.join(config.plotsbasedir, "reweightqqZZ", str(analysis), str(category))
  mkdir_p(saveasdir)
  for ext in "png eps root pdf".split():
    c.SaveAs(os.path.join(saveasdir, "{}.{}".format(disc, ext)))

if __name__ == "__main__":
  for analysis in analyses:
    for category in categories:
      discs = TemplatesFile("2e2mu", "ggh", category, analysis).discriminants
      for disc in discs:
        plot(disc.name, category, analysis)
