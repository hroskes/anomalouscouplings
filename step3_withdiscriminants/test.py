#!/usr/bin/env python

assert __name__ == "__main__"

from helperstuff import config
from helperstuff import style
from helperstuff.enums import *
from helperstuff.plotfromtree import plotfromtree
from helperstuff.samples import *
import ROOT
import os

#========================
#inputs
#weight, bins, min, max, category can be None
productionmode = "WH"
disc           = "D_L1Zgint_VBF"
weight         = ReweightingSample(productionmode, "fL1proddec0.5")
bins           = None
min            = None
max            = None

enrich         = False
masscut        = True
normalizeto1   = False

channel        = "2e2mu"
category       = None

skip           = []
#========================

hstack = ROOT.THStack()
legend = ROOT.TLegend(.6, .5, .9, .9)
cache = []

ROOT.gErrorIgnoreLevel = ROOT.kError

c = ROOT.TCanvas()
hs = {}

productionmode = ProductionMode(productionmode)

def hypothesestouse():
    for hypothesis in hypotheses:
        if hypothesis not in productionmode.generatedhypotheses: continue
        if productionmode == "ggH" and hypothesis == "L1": continue
        if skip is not None and hypothesis in skip: continue
        yield hypothesis

for hypothesis in hypotheses:
    if hypothesis not in hypothesestouse():
        try:
            os.remove(os.path.join(config.plotsbasedir, "TEST", "reweighting", "{}.png".format(hypothesis)))
        except OSError:
            pass

for color, hypothesis in enumerate(hypothesestouse(), start=1):
    if color == 5: color = ROOT.kYellow+3

    hname = "h{}".format(hypothesis)

    h = hs[hypothesis] = plotfromtree(
      productionmode=productionmode,
      hypothesis=hypothesis,
      weight=weight,
      disc=disc,
      bins=bins,
      min=min,
      max=max,
      enrich=enrich,
      masscut=masscut,
      normalizeto1=normalizeto1,
      channel=channel,
      category=category,
      color=color,
      hname=hname,
    )

    hstack.Add(h)
    cache.append(h)
    legend.AddEntry(h, str(hypothesis), "l")
    print "{:10} {:8.3g}      {:8.3g}".format(hypothesis, h.Integral(), h.GetEffectiveEntries())
    try:
      os.makedirs(os.path.join(config.plotsbasedir, "TEST", "reweighting"))
    except OSError:
      pass
    for ext in "png eps root pdf".split():
      c.SaveAs(os.path.join(config.plotsbasedir, "TEST", "reweighting", "{}.{}".format(hypothesis, ext)))

hstack.Draw("hist")
hstack.GetXaxis().SetTitle(h.GetXaxis().GetTitle())
for ext in "png eps root pdf".split():
  c.SaveAs(os.path.join(config.plotsbasedir, "TEST", "reweighting", "test.{}".format(ext)))

