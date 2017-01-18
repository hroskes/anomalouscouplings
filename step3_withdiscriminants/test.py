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
productionmode = "VBF"
disc           = "D_CP_VBF"
weight         = ArbitraryCouplingsSample('VBF', g1=1, g2=0, g4=-0.297979, g1prime2=0)
bins           = None
min            = None
max            = None

enrich         = False
masscut        = True
normalizeto1   = True

channel        = "2e2mu"
category       = None
#========================

hstack = ROOT.THStack()
legend = ROOT.TLegend(.6, .5, .9, .9)
cache = []

ROOT.gErrorIgnoreLevel = ROOT.kError

c = ROOT.TCanvas()
hs = {}

def hypothesestouse():
    for hypothesis in hypotheses:
        if productionmode == "ggH" and hypothesis not in decayonlyhypotheses: continue
        if productionmode in ["VBF", "ZH", "WH"] and hypothesis not in prodonlyhypotheses: continue
        if productionmode in ["HJJ", "ttH"] and hypothesis not in hffhypotheses: continue
        yield hypothesis

for hypothesis in hypotheses:
    if hypothesis not in hypothesestouse():
        try:
            os.remove(os.path.join(config.plotsbasedir, "TEST", "{}.png".format(hypothesis)))
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
    print "{:10} {:.3g}".format(hypothesis, h.Integral())
    try:
      os.makedirs(os.path.join(config.plotsbasedir, "TEST"))
    except OSError:
      pass
    for ext in "png eps root pdf".split():
      c.SaveAs(os.path.join(config.plotsbasedir, "TEST", "{}.{}".format(hypothesis, ext)))

hstack.Draw("histnostack")
hstack.GetXaxis().SetTitle(h.GetXaxis().GetTitle())
for ext in "png eps root pdf".split():
  c.SaveAs(os.path.join(config.plotsbasedir, "TEST", "test.{}".format(ext)))

"""
hint = hs["0+"].Clone("hint")
hint.Add(hs["0-"])
hint.Scale(-.5)
hint.Add(hs["fa3prod0.5"])

hint.Draw("hist")
#c.SaveAs(os.path.join(config.plotsbasedir, "TEST", "a1a3int.png"))
"""
