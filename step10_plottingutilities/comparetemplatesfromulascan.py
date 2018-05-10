#!/usr/bin/env python

import argparse
import os

import ROOT

from helperstuff import config
from helperstuff.copyplots import copyplots
from helperstuff.enums import Analysis, categories, channels, ProductionMode
from helperstuff.samples import Sample
from helperstuff.plotfromtree import plotfromtree
from helperstuff.templates import Template
from helperstuff.utilities import mkdir_p

def compare(template, axis):
  c = ROOT.TCanvas()

  t = template.gettemplate()
  h = (t.ProjectionY, t.ProjectionZ, t.ProjectionX)[axis]()
  h.Scale(1/h.Integral())
  h.SetLineColor(2)

  flavor = str(template.channel) if template.productionmode in ("ggZZ", "VBF bkg") else None

  s = Sample(template.reweightingsampleplus, template.production, flavor)
  h2 = plotfromtree(
    reweightfrom=s,
    disc=template.discriminants[axis],
    normalizeto1=True,
    channel=template.channel,
    category=template.category,
    analysis=template.analysis,
    color=4,
  )

  hstack = ROOT.THStack()
  hstack.Add(h)
  hstack.Add(h2)
  hstack.Draw("hist nostack")
  hstack.GetXaxis().SetTitle(h2.GetXaxis().GetTitle())

  saveas = os.path.join(config.plotsbasedir, "templateprojections", "comparetoUlascan", str(template.analysis), str(template.category), str(template.channel), template.discriminants[axis].name+"_"+str(template.productionmode)+"_"+str(template.hypothesis))
  saveas = saveas.replace("_None", "")
  mkdir_p(os.path.dirname(saveas))
  for ext in "png eps root pdf".split():
    c.SaveAs(saveas+"."+ext)

if __name__ == "__main__":
  try:
    for channel in channels:
      for analysis in Analysis.items(lambda x: x in ("fa3", "fa2", "fL1")):
        for category in categories:
          for productionmode in ProductionMode.items(lambda x: x in ("ggH", "VBF", "ZH", "WH", "ggZZ", "qqZZ")):
            if productionmode.issignal: hypotheses = analysis.purehypotheses
            else: hypotheses = None,
            for hypothesis in hypotheses:
              for axis in 0, 1, 2:
                compare(Template(analysis, productionmode, hypothesis, category, channel, "180416"), axis)


  finally:
    copyplots("templateprojections/comparetoUlascan")
