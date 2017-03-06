#!/usr/bin/env python

from itertools import product
import os

import ROOT

from helperstuff import config, style
from helperstuff.discriminants import discriminant
from helperstuff.enums import analyses, categories, Category, hffhypotheses
from helperstuff.samples import ReweightingSample, ReweightingSamplePlus, Sample
from helperstuff.plotfromtree import Line, plotfromtree
from helperstuff.templates import TemplatesFile
from helperstuff.utilities import tfiles

assert len(config.productionsforcombine) == 1
production = config.productionsforcombine[0]
color = 0
def makeline(sample, title, reweightfrom=None):
  global color
  color += 1
  thecolor = color
  assert thecolor < 8
  if thecolor == 7: thecolor = ROOT.kViolet-1
  if thecolor == 5: thecolor = 6
  return Line(sample, title, thecolor, reweightfrom)

def lines(productionmode, analysis):
  global color
  color = 0
  if productionmode == "HJJ":
    yield makeline(ReweightingSample(productionmode, "0+", "Hff0+"), "ggH SM JHUGen H+jj")
    yield makeline(ReweightingSample(productionmode, "0+", "Hff0-"), "ggH 0^{-} JHUGen H+jj")
    yield makeline(ReweightingSample(productionmode, "0+", "fCP0.5"), "ggH f_{CP}=0.5 JHUGen H+jj")
    yield makeline(ReweightingSamplePlus("ggH", "0+", "MINLO"), "ggH SM MINLO")
    yield makeline(ReweightingSamplePlus("ggH", "0+", "NNLOPS"), "ggH SM NNLOPS")
  elif productionmode == "ttH":
    yield makeline(ReweightingSample(productionmode, "0+", "Hff0+"), "ttH SM JHUGen")
    yield makeline(ReweightingSample(productionmode, "0+", "Hff0-"), "ttH 0^{-} JHUGen")
    yield makeline(ReweightingSample(productionmode, "0+", "fCP0.5"), "ttH f_{CP}=0.5 JHUGen")
  else:
    reweightfrom = None
    if analysis == "fL1Zg": reweightfrom = ReweightingSample(productionmode, "0+")

    yield makeline(ReweightingSample(productionmode, analysis.purehypotheses[0]), "{} SM JHUGen".format(productionmode), reweightfrom)

    if analysis == "fL1Zg":
      reweightfrom = ReweightingSample(productionmode, "L1")
    if analysis == "fL1Zg" and productionmode == "WH":
      color += 2
    else:
      yield makeline(ReweightingSample(productionmode, analysis.purehypotheses[1]), "{} {}=1 JHUGen".format(productionmode, analysis.title(superscript=productionmode)), reweightfrom)
      yield makeline(ReweightingSample(productionmode, analysis.mixprodhypothesis), "{} {}=0.5 JHUGen".format(productionmode, analysis.title(superscript=productionmode)), reweightfrom)

  if productionmode == "HJJ":
    yield makeline(ReweightingSample("ggH", "SM"), "ggH POWHEG")
  elif productionmode == "WH":
    yield makeline(Sample("WplusH", "SM", "POWHEG", production), "W^{+}H SM POWHEG")
    yield makeline(Sample("WminusH", "SM", "POWHEG", production), "W^{-}H SM POWHEG")
  elif productionmode == "ttH":
    yield makeline(Sample(productionmode, "SM", "Hff0+", "POWHEG", production), "ttH SM POWHEG")
  else:
    yield makeline(Sample(productionmode, "SM", "POWHEG", production), "{} SM POWHEG".format(productionmode))

cache = {}

for analysis, category in product(analyses, categories):
  if category == "Untagged": continue
  for discname, title, bins, min, max in TemplatesFile("ggh", "2e2mu", analysis, category, production).discriminants[0:2]:
    for productionmode in "HJJ", "ttH", "VBF", "ZH", "WH":
      categoryname = "category_"+analysis.categoryname
      assert len(category.idnumbers) == 1
      categoryid = category.idnumbers[0]
      hstack = ROOT.THStack(productionmode, title)
      legend = ROOT.TLegend(.6, .7, .9, .9)
      legend.SetBorderSize(0)
      legend.SetFillStyle(0)
      for sample, title, color, reweightfrom in lines(productionmode, analysis):
        hffhypothesis = None
        if productionmode in ("HJJ", "ttH"): hffhypothesis = reweightfrom.hffhypothesis
        h = cache[analysis,discname,sample] = plotfromtree(disc=discname, color=color, reweightto=sample, reweightfrom=reweightfrom, normalizeto1=True, category=category, analysis=analysis)
        hstack.Add(h)
        legend.AddEntry(h, title, "l")

      c = ROOT.TCanvas()
      hstack.Draw("hist nostack")
      hstack.GetXaxis().SetTitle(h.GetXaxis().GetTitle())
      legend.Draw()

      directory = os.path.join(config.plotsbasedir, "templateprojections", "compare_POWHEG_JHUGen", str(analysis), str(productionmode))
      try:
        os.makedirs(directory)
      except OSError:
        pass
      for ext in "png eps root pdf".split():
        c.SaveAs(os.path.join(directory, "{}.{}".format(discname, ext)))
