#!/usr/bin/env python

import os

import ROOT

from helperstuff import config, style
from helperstuff.combinehelpers import gettemplate
from helperstuff.utilities import mkdir_p, PlotCopier

def draw(channel, production, plotcopier=ROOT):
  store = []
  hstack = ROOT.THStack("D_bkg", "D_bkg")
  l = ROOT.TLegend(.6, .6, .9, .9)
  l.SetBorderSize(0)
  l.SetFillStyle(0)
  for color, syst in enumerate(("", "ScaleUp", "ScaleDown", "ResUp", "ResDown"), start=1):
    title = syst if syst else "Nominal"
    if color == 5: color = 6
    template = gettemplate("ggH", "fa3", channel, "Untagged", production, "0+", syst)
    h = template.ProjectionZ(template.GetName()+syst)
    h.SetLineColor(color)
    hstack.Add(h)
    store.append(h)
    l.AddEntry(h, title, "l")
  c = plotcopier.TCanvas()
  hstack.Draw("hist nostack")
  hstack.GetXaxis().SetTitle("D_{bkg}")
  l.Draw()
  saveasdir = os.path.join(config.plotsbasedir, "xchecks", "m4lsystematics", str(production.year))
  mkdir_p(saveasdir)
  for ext in "png eps root pdf".split():
    c.SaveAs(os.path.join(saveasdir, "{}.{}".format(channel, ext)))

if __name__ == "__main__":
  with PlotCopier() as pc:
    for c in "2e2mu", "4e", "4mu":
      for production in config.productionsforcombine:
        draw(c, production, pc)
