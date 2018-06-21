#!/usr/bin/env python

import argparse

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--yields", action="store_true")
  args = parser.parse_args()

import os

import ROOT

from helperstuff import config, style
from helperstuff.copyplots import copyplots
from helperstuff.enums import analyses, Analysis, categories, channels, ProductionMode, productions
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

def compareyields(analysis, productionmode, production):
  from helperstuff.combinehelpers import getrate
  productionmode = ProductionMode(productionmode)
  hypothesis = "0+" if productionmode.issignal else None
  for category in categories:
    for channel in channels:
      r = Template(analysis, productionmode, category, channel, str(production)+"_Ulascan", hypothesis).gettemplate().Integral()
      r2 = getrate(analysis, productionmode, category, channel, production, 1)
      print "{:10} {:5} {:8.4f} {:8.4f} {:.2%}".format(category, channel, r, r2, r/r2-1)

if __name__ == "__main__":
  if args.yields:
    for analysis in analyses:
      for productionmode in ProductionMode.items(lambda x: x in ("ggH", "VBF", "ZH", "WH", "ggZZ", "qqZZ")):
        for production in productions:
          if analysis != "fa3": continue
          print "================="
          print "{:5} {:>4} {:6}".format(analysis, productionmode, production)
          print "================="
          compareyields(analysis, productionmode, production)
  else:
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
