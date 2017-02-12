#!/usr/bin/env python

import os

import ROOT

from helperstuff import config, style
from helperstuff.discriminants import discriminant
from helperstuff.enums import analyses, Analysis
from helperstuff.samples import ArbitraryCouplingsSample, ReweightingSample, Sample
from helperstuff.plotfromtree import Line, plotfromtree
from helperstuff.rootoverloads import histogramadd

assert len(config.productionsforcombine) == 1
production = config.productionsforcombine[0]

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
      color=color*2-2+linestyle
      if cut is None: continue
#      if cut is not None: continue

      htotal = hs[disc,sample,style] = None
      for hypothesis in sample.productionmode.generatedhypotheses:
        reweightfrom = Sample(hypothesis, sample.productionmode, production)
        h = hs[disc,sample,style,reweightfrom] = plotfromtree(
          reweightfrom=reweightfrom,
          reweightto=sample,
          disc=disc,
          normalizeto1=False,
          color=color,
          cut=cut,
#          linestyle=linestyle,
          channel=channel,
          category="VBFtagged",
          analysis="fa2",
          bins=20,
        )
        if htotal is None: htotal = h
        else: htotal += h
      htotal.Scale(1/htotal.Integral())
      htotal.SetLineWidth(2)
      htotal.SetMarkerStyle(20+linestyle)
      hstack.Add(htotal)
      usetitle = title
      if cut: usetitle += " ("+cut.replace("ZZP", "p")+")"
      l.AddEntry(htotal, usetitle, "l")

  hstack.Draw("PEnostack")
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

