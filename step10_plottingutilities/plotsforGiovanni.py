#!/usr/bin/env python

import os

import ROOT

from helperstuff import config, style
from helperstuff.samples import ReweightingSample
from helperstuff.discriminants import discriminant
from helperstuff.plotfromtree import plotfromtree
from helperstuff.utilities import mkdir_p, PlotCopier

def hist(hypothesis, disc):
  return plotfromtree(
    reweightfrom=ReweightingSample("ggH", "0+"),
    reweightto=ReweightingSample("ggH", hypothesis.hypothesis),
    production="190821_2016",
    disc=disc,
    normalizeto1=True,
    color=hypothesis.ffHlinecolor,
    linestyle=hypothesis.ffHlinestyle,
  )

from discriminantplots import HypothesisLine, purehypothesislines, a3mixdecay, a2mixdecay, L1mixdecay
bestfit = HypothesisLine("BestFit19009", 1, 1, 1, 2, "best fit")

def plot(disc, pc=ROOT):
  print disc
  disc = discriminant(disc)

  hypothesislines = purehypothesislines[:]
  if "0hplus" in disc.name or "int" in disc.name: hypothesislines.append(a2mixdecay)
  if "CP" in disc.name: hypothesislines.append(a3mixdecay)
  if "L1" in disc.name: hypothesislines.append(L1mixdecay)
  hypothesislines.append(bestfit)

  hists = [hist(h, disc) for h in hypothesislines]

  hstack = ROOT.THStack(disc.name, "")
  l = ROOT.TLegend(.6, .6, .9, .9)
  l.SetFillStyle(0)
  l.SetBorderSize(0)
  for hypo, h in zip(hypothesislines, hists):
    hstack.Add(h)
    l.AddEntry(h, hypo.legendname, "l")
  c = pc.TCanvas()
  hstack.Draw("hist nostack")
  hstack.GetXaxis().SetTitle(disc.title)
  l.Draw()
  mkdir_p(os.path.join(config.plotsbasedir, "templateprojections", "bestfitcheck"))
  for ext in "png pdf root C".split():
    c.SaveAs(os.path.join(config.plotsbasedir, "templateprojections", "bestfitcheck", disc.name+"."+ext))

if __name__ == "__main__":
  with PlotCopier() as pc:
    plot("D_0minus_decay_3bins", pc)
    plot("D_0hplus_decay_3bins", pc)
    plot("D_L1_decay_3bins", pc)
    plot("D_L1Zg_decay_3bins", pc)
    plot("D_CP_decay_2bins", pc)
    plot("D_int_decay_2bins", pc)
    plot("D_bkg_3bins", pc)
