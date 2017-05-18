#!/usr/bin/env python

import os

import ROOT

from helperstuff import config, style
from helperstuff.combinehelpers import gettemplate
from helperstuff.utilities import mkdir_p

assert len(config.productionsforcombine) == 1
production = config.productionsforcombine[0]

store = []

def draw(channel):
  hstack = ROOT.THStack("D_bkg", "D_bkg")
  l = ROOT.TLegend(.6, .6, .9, .9)
  l.SetBorderSize(0)
  l.SetFillStyle(0)
  for color, syst in enumerate(("", "ScaleUp", "ScaleDown", "ResUp", "ResDown"), start=1):
#  for color, syst in enumerate(("", "MINLOUp", "MINLODn"), start=1):
    title = syst if syst else "Nominal"
    if color == 5: color = 6
    h = gettemplate("ggH", "fa3", channel, "Untagged", production, "0+", syst).ProjectionZ()
    h.SetLineColor(color)
    hstack.Add(h)
    store.append(h)
    l.AddEntry(h, title, "l")
  c = ROOT.TCanvas()
  hstack.Draw("hist nostack")
  hstack.GetXaxis().SetTitle("D_{bkg}")
  l.Draw()
  saveasdir = os.path.join(config.plotsbasedir, "xchecks", "m4lsystematics")
  mkdir_p(saveasdir)
  for ext in "png eps root pdf".split():
    c.SaveAs(os.path.join(saveasdir, "{}.{}".format(channel, ext)))

if __name__ == "__main__":
  for c in "2e2mu", "4e", "4mu":
    draw(c)
