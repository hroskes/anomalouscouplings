#!/usr/bin/env python

import os

import ROOT

from helperstuff import config, style
from helperstuff.discriminants import discriminant
from helperstuff.enums import analyses, Analysis
from helperstuff.samples import ArbitraryCouplingsSample, ReweightingSample
from helperstuff.plotfromtree import Line, plotfromtree


def makeplot(disc, channel):
  lines = (
           Line(ReweightingSample("ggH", "fa2-0.5"), "ggH", 1, ReweightingSample("ggH", "a2")),
           Line(ReweightingSample("VBF", "fa2-0.5"), "VBF", 2, ReweightingSample("VBF", "a2")),
#           Line(ReweightingSample("ZH", "fa2-0.5"), "ZH", 4, ReweightingSample("ZH", "a2")),
#           Line(ReweightingSample("WH", "fa2-0.5"), "WH", ROOT.kGreen+3, ReweightingSample("WH", "a2")),
          )

  hs = {}
  hstack = ROOT.THStack(disc, "")
  l = ROOT.TLegend(.15, .5, .4, .9)
  l.SetBorderSize(0)
  l.SetFillStyle(0)
  c = ROOT.TCanvas()

  for sample, title, color, reweightfrom in lines:
    for linestyle, cut in enumerate((None, "ZZPt<200", "ZZPt>200"), start=1):
#      color=linestyle
#      if cut is None: continue
      if cut is not None: continue
      linestyle=1
#      if cut is not None: continue
      h = hs[disc,sample,style] = plotfromtree(
        reweightfrom=reweightfrom,
        reweightto=sample,
        disc=disc,
        normalizeto1=True,
        color=color,
        cut=cut,
        linestyle=linestyle,
        channel=channel,
        category="VBFtagged",
        analysis="fa2",
        bins=20,
      )
      h.SetLineWidth(2)
      hstack.Add(h)
      usetitle = title
      if cut: usetitle += " ("+cut.replace("ZZP", "p")+")"
      l.AddEntry(h, usetitle, "l")

  hstack.Draw("hist nostack")
  l.Draw()

  saveasdir = os.path.join(config.plotsbasedir, "fa2_pTeffect", str(channel))
  try:
    os.makedirs(saveasdir)
  except OSError:
    pass
  for ext in "png eps root pdf".split():
    c.SaveAs(os.path.join(saveasdir, "{}.{}".format(disc, ext)))

if __name__ == "__main__":
  for channel in "2e2mu", "4e", "4mu":
    makeplot("D_0hplus_decay", channel)
    makeplot("D_int_decay", channel)

