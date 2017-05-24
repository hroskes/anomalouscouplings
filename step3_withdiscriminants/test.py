#!/usr/bin/env python

assert __name__ == "__main__"

import itertools
import os

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
normalizeto1 = False
color = 1

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
      cut=cut
    )
    h.SetMinimum(0)

    cache.append(h)
    print "{:20} {:20} {:8.3g}      {:8.3g}".format(disc, disc2, h.Integral(), h.GetEffectiveEntries())
    try:
      os.makedirs(os.path.join(config.plotsbasedir, "TEST", "reweighting"))
    except OSError:
      pass
    for ext in "png eps root pdf".split():
      c.SaveAs(os.path.join(config.plotsbasedir, "TEST", "{}{}.{}".format(disc, disc2, ext)))

copyplots("TEST")
