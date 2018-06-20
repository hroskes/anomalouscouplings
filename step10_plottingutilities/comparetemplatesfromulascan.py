#!/usr/bin/env python

import argparse
import os

import ROOT

from helperstuff import config, style
from helperstuff.copyplots import copyplots
from helperstuff.enums import Analysis, categories, channels, ProductionMode
from helperstuff.samples import Sample
from helperstuff.plotfromtree import plotfromtree
from helperstuff.templates import Template
from helperstuff.utilities import mkdir_p

c = ROOT.TCanvas()

def compare(template, axis):
  t = Template(str(template.production)+"_Ulascan", template.productionmode, template.hypothesis, template.category, template.channel, template.analysis, template.shapesystematic).gettemplate()
  h = (t.ProjectionY, t.ProjectionZ, t.ProjectionX)[axis]()
  h.Scale(1/h.Integral())
  h.SetLineColor(2)

  t2 = template.gettemplate()
  h2 = (t2.ProjectionX, t2.ProjectionY, t2.ProjectionZ)[axis]()
  h2.Scale(1/h2.Integral())
  h2.SetLineColor(4)

  nbinstomerge = 1.0 * h.GetNbinsX() / h2.GetNbinsX()
  assert nbinstomerge == int(nbinstomerge)
  if nbinstomerge != 1:
    nbinstomerge = int(nbinstomerge)
    h.Rebin(nbinstomerge)

  hstack = ROOT.THStack()
  hstack.Add(h)
  hstack.Add(h2)
  hstack.Draw("hist nostack")
  hstack.GetXaxis().SetTitle(template.discriminants[axis].title)

  l = ROOT.TLegend(.6, .7, .9, .9)
  l.SetBorderSize(0)
  l.SetFillStyle(0)
  l.AddEntry(h, "Ulascan", "l")
  l.AddEntry(h2, "Heshy", "l")
  l.Draw()

  saveas = os.path.join(config.plotsbasedir, "templateprojections", "comparetoUlascan", str(template.production.year), str(template.analysis), str(template.category), str(template.channel), template.discriminants[axis].name+"_"+str(template.productionmode)+"_"+str(template.hypothesis))
  saveas = saveas.replace("_None", "")
  mkdir_p(os.path.dirname(saveas))
  for ext in "png eps root pdf".split():
    c.SaveAs(saveas+"."+ext)

if __name__ == "__main__":
  try:
    for channel in channels:
      for analysis in Analysis.items(lambda x: x in ("fa3", "fa2", "fL1", "fL1Zg")):
        for category in categories:
          for productionmode in ProductionMode.items(lambda x: x in ("ggH", "VBF", "ZH", "WH", "ggZZ", "qqZZ")):
            if productionmode == "WH" and analysis == "fL1Zg": continue
            if productionmode.issignal: hypotheses = analysis.purehypotheses
            elif productionmode == "qqZZ": hypotheses = "ext",
            else: hypotheses = None,
            for hypothesis in hypotheses:
              for production in "180530", "180531":
                for axis in 0, 1, 2:
                  compare(Template(analysis, productionmode, hypothesis, category, channel, production), axis)


  finally:
    copyplots("templateprojections/comparetoUlascan")
