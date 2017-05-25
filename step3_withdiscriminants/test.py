#!/usr/bin/env python

assert __name__ == "__main__"

import itertools
import math
import os
import random

import ROOT

from helperstuff import config
from helperstuff import style
from helperstuff.copyplots import copyplots
from helperstuff.enums import *
from helperstuff.plotfromtree import plotfromtree
from helperstuff.samples import *

cache = []

ROOT.gErrorIgnoreLevel = ROOT.kError

c = ROOT.TCanvas()
hs = {}

productionmode = "ggH"
hypothesis = "0+_photoncut"
reweightto = None
bins = min = max = bins2 = min2 = max2 = channel = category = analysis = hname = cut = None
enrich = False
masscut = True
normalizeto1 = True
color = 1

def mutualinformation(h2):
  h2 = h2.Clone("h{}".format(random.random()))
  h2.Scale(1/h2.Integral())
  hx = h2.ProjectionX()
  hy = h2.ProjectionY()

  MI = 0
  #https://en.wikipedia.org/wiki/Mutual_information#Definition
  for x in range(1, h2.GetNbinsX()+1):
    for y in range(1, h2.GetNbinsY()+1):
      MI += h2.GetBinContent(x, y) * math.log(h2.GetBinContent(x, y) / (hx.GetBinContent(x) * hy.GetBinContent(y)))
  return MI

discriminants = "D_L1_decay", "D_L1int_decay", "D_L1Zg_decay", "D_L1Zgint_decay", "D_L1L1Zg_decay", "D_L1L1Zgint_decay"

for (i1, disc), (i2, disc2) in itertools.product(enumerate(discriminants), enumerate(discriminants)):
    if i1 >= i2: continue

    hname = "h{}{}".format(disc, disc2)

    h = hs[hypothesis] = plotfromtree(
      reweightfrom=ReweightingSample(productionmode, hypothesis),
      reweightto=reweightto,
      disc=disc,
      bins=bins,
      min=min,
      max=max,
      disc2=disc2,
      bins2=bins2,
      min2=min2,
      max2=max2,
      enrich=enrich,
      masscut=masscut,
      normalizeto1=normalizeto1,
      channel=channel,
      category=category,
      analysis=analysis,
      color=color,
      hname=hname,
      cut=cut,
    )
    h.SetMinimum(0)

    cache.append(h)
    print "{:20} {:20}      {:8.2f}".format(disc, disc2, mutualinformation(h))
    try:
      os.makedirs(os.path.join(config.plotsbasedir, "TEST", "reweighting"))
    except OSError:
      pass
    for ext in "png eps root pdf".split():
      c.SaveAs(os.path.join(config.plotsbasedir, "TEST", "{}{}.{}".format(disc, disc2, ext)))

for disc in discriminants:
  hstack = ROOT.THStack("hstack{}".format(disc), "")
  l = ROOT.TLegend(.6, .7, .9, .9)
  for color, (hypothesis, title) in enumerate(zip(("0+_photoncut", "fL10.5_photoncut", "L1Zg"), ("SM", "f_{#Lambda1}=0.5", "#Lambda_{1}^{Z#gamma}")), start=1):
    hname = "h{}{}".format(disc, hypothesis)

    h = plotfromtree(
      reweightfrom=ReweightingSample(productionmode, hypothesis),
      reweightto=reweightto,
      disc=disc,
      bins=bins,
      min=min,
      max=max,
      enrich=enrich,
      masscut=masscut,
      normalizeto1=normalizeto1,
      channel=channel,
      category=category,
      analysis=analysis,
      color=color,
      hname=hname,
      cut=cut,
    )
    h.SetMinimum(0)
    l.AddEntry(h, title, "l")

    cache.append(h)
    hstack.Add(h)

  hstack.Draw("hist nostack")
  l.Draw()
  try:
    os.makedirs(os.path.join(config.plotsbasedir, "TEST", "reweighting"))
  except OSError:
    pass
  for ext in "png eps root pdf".split():
    c.SaveAs(os.path.join(config.plotsbasedir, "TEST", "{}.{}".format(disc, ext)))

copyplots("TEST")
