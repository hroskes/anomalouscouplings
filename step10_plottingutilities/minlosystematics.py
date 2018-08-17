#!/usr/bin/env python

import os

import ROOT

from helperstuff import config, style
from helperstuff.combinehelpers import gettemplate
from helperstuff.templates import Template
from helperstuff.utilities import mkdir_p, PlotCopier

def draw(channel, category, production, plotcopier=ROOT):
  store = []
  for i, axis in enumerate((ROOT.TH3.ProjectionX, ROOT.TH3.ProjectionY, ROOT.TH3.ProjectionZ)):
    hstack = ROOT.THStack("D_bkg", "D_bkg")
    l = ROOT.TLegend(.6, .6, .9, .9)
    l.SetBorderSize(0)
    l.SetFillStyle(0)
    for color, syst in enumerate(("", "MINLOUp", "MINLODown"), start=1):
      title = syst if syst else "Nominal"
      if color == 5: color = 6
      template = Template("ggH", analysis, channel, category, production, "0+", syst)
      disc = template.discriminants[i]
      h = axis(template.gettemplate(), template.templatename()+disc.name+syst)
      h.SetLineColor(color)
      hstack.Add(h)
      store.append(h)
      l.AddEntry(h, title, "l")
    c = plotcopier.TCanvas()
    hstack.Draw("hist nostack")
    hstack.GetXaxis().SetTitle(disc.title)
    l.Draw()
    saveasdir = os.path.join(config.plotsbasedir, "xchecks", "minlosystematics", str(production.year), str(analysis), str(channel))
    mkdir_p(saveasdir)
    for ext in "png eps root pdf".split():
      c.SaveAs(os.path.join(saveasdir, "{}.{}".format(disc.name, ext)))

if __name__ == "__main__":
  with PlotCopier() as pc:
    for ch in "2e2mu", "4e", "4mu":
      for ca in "VBFtagged", "VHHadrtagged":
        for analysis in "fa3", "fa2", "fL1", "fL1Zg":
          for production in config.productionsforcombine:
            draw(ch, ca, production, plotcopier=pc)
