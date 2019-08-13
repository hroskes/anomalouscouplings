#!/usr/bin/env python

import os

import ROOT

from helperstuff import config, style
from helperstuff.combinehelpers import gettemplate
from helperstuff.plotfromtree import plotfromtree
from helperstuff.samples import ReweightingSample
from helperstuff.utilities import mkdir_p, PlotCopier

def draw(channel, category, production, plotcopier=ROOT):
  store = []
  hstack = ROOT.THStack("D_bkg", "D_bkg")
  l = ROOT.TLegend(.6, .6, .9, .9)
  l.SetBorderSize(0)
  l.SetFillStyle(0)

  for color, syst in enumerate(("", "ScaleUp", "ScaleDown", "ResUp", "ResDown"), start=1):
    title = syst if syst else "Nominal"
    if color == 5: color = 6

    disc = "D_bkg_3bins"
    if category == "VBFtagged": disc ="D_bkg_VBFdecay_3bins"
    if category == "VHHadrtagged": disc ="D_bkg_HadVHdecay_3bins"
    if syst: disc = disc.replace("_3bins", "_"+syst+"_3bins")

    h = plotfromtree(
      reweightfrom=ReweightingSample("ggH", "0+"),
      production=production,
      channel=channel,
      disc=disc,
      enrich=False,
      category=category,
      masscut=True,
      analysis="fa3fa2fL1fL1Zg_morecategories",
    )

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
  for ext in "png C root pdf".split():
    c.SaveAs(os.path.join(saveasdir, "{}{}.{}".format(category, channel, ext)))

if __name__ == "__main__":
  with PlotCopier() as pc:
    for c in "2e2mu", "4e", "4mu":
      for ca in "Untagged", "VBFtagged", "VHHadrtagged", "VHLepttagged", "VBF1jtagged", "Boosted":
        for production in config.productionsforcombine:
          draw(channel=c, category=ca, production=production, plotcopier=pc)
