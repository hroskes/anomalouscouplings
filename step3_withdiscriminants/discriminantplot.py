#!/usr/bin/env python

assert __name__ == "__main__"

from helperstuff import config
from helperstuff import style
from helperstuff.CJLSTscripts import *
from helperstuff.enums import *
from helperstuff.plotfromtree import plotfromtree
from helperstuff.samples import *
from helperstuff.utilities import PlotCopier
import ROOT
import os

#========================
#inputs
productionmode = "VBF"
disc           = "ZZPt"
reweightto     = None
bins           = 3
min            = 40
max            = 280

production     = "190821_2018"

enrich         = False
masscut        = True
normalizeto1   = True

channel        = None

category       = "Untagged"
analysis       = "fa3fa2fL1fL1Zg"

cut            = None

disc2          = None
bins2          = None
min2           = None
max2           = None

skip           = []

logscale       = False
setminimum     = 1e-7
#========================

hstack = ROOT.THStack()
legend = ROOT.TLegend(.6, .5, .9, .9)
cache = []

ROOT.gErrorIgnoreLevel = ROOT.kError

productionmode = ProductionMode(productionmode)

def hypothesestouse():
  if reweightto is not None:
    for hypothesis in productionmode.generatedhypotheses(production):
      yield hypothesis
    return
  if productionmode == "VBF":
    for hypothesis in hypotheses:
      if not hypothesis.ispure: continue
      if hypothesis in ("0+", "0-", "a2", "L1", "L1Zg", "fa3prod0.5", "fa2prod0.5", "fL1Zgprod0.5", "fL1prod0.5"): yield hypothesis
  elif productionmode == "ZH":
    for hypothesis in hypotheses:
      if hypothesis in ("0+", "0-", "a2", "L1", "L1Zg", "fa3prod0.5", "fa2prod0.5", "fL1Zgprod0.5", "fL1prod0.5"): yield hypothesis
  else:
    for hypothesis in hypotheses:
      if hypothesis in ("0+", "0-", "a2", "L1", "L1Zg", "fa30.5", "fa2-0.5", "fL1Zg0.5", "fL10.5"): yield hypothesis

def hffhypothesistouse():
  if productionmode == "ttH":
    return "Hff0+"
  else:
    return None

with PlotCopier() as pc:
  c = pc.TCanvas()
  if logscale: c.SetLogy()
  hs = {}

  for hypothesis in hypotheses:
    if hypothesis not in list(hypothesestouse()):
      try:
        pc.remove(os.path.join(config.plotsbasedir, "TEST", "reweighting", "{}.png".format(hypothesis)))
      except OSError:
        pass

  for color, hypothesis in enumerate(hypothesestouse(), start=1):
    if color == 5: color = ROOT.kYellow+3

    hname = "h{}".format(hypothesis)
    rwtfrom = hypothesis

    h = hs[hypothesis] = plotfromtree(
      reweightfrom=ReweightingSample(productionmode, rwtfrom, hffhypothesistouse()),
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
      production=production if hypothesis != "0+" else "190821_2016",
    )
    assert production == "190821_2018", "^^^^^^^^^^^"
    if not logscale: h.SetMinimum(0)
    if setminimum is not None: h.SetMinimum(setminimum)

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

  if reweightto is None:
    option = "hist nostack"
  else:
    option = "hist"
  hstack.Draw(option)
  hstack.GetXaxis().SetTitle(h.GetXaxis().GetTitle())
  for ext in "png eps root pdf".split():
    c.SaveAs(os.path.join(config.plotsbasedir, "TEST", "reweighting", "test.{}".format(ext)))
