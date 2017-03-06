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
productionmode = "ttH"
disc           = "D_CP_VBF"
reweightto     = ReweightingSample(productionmode, "0+", "Hff0+")
bins           = None
min            = None
max            = None

enrich         = False
masscut        = True
normalizeto1   = False

channel        = None

category       = "VBFtagged"
analysis       = "fa3"

cut            = "D_bkg<1./40"

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
        if skip is not None and hypothesis in skip: continue
        yield hypothesis

def hffhypothesistouse():
    if productionmode == "ttH":
        return "Hff0+"
    else:
        return None

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
      reweightfrom=ReweightingSample(productionmode, hypothesis, hffhypothesistouse()),
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
      cut=cut
    )
    h.SetMinimum(0)

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

