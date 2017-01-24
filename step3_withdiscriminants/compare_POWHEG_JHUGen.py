#!/usr/bin/env python

from itertools import product
import os

import ROOT

from helperstuff import config, style
from helperstuff.discriminants import discriminant
from helperstuff.enums import analyses, categories, Category, hffhypotheses
from helperstuff.samples import Sample
from helperstuff.templates import TemplatesFile
from helperstuff.utilities import tfiles

assert len(config.productionsforcombine) == 1
production = config.productionsforcombine[0]

def samples(productionmode):
  if productionmode in ("HJJ", "ttH"):
    for hffhypothesis in hffhypotheses:
      yield Sample(productionmode, hffhypothesis, production, "0+")
  else:
    for hypothesis in ("0+", "0-", "fa3prod0.5"):
      yield Sample(productionmode, hypothesis, production)
  if productionmode == "HJJ":
    yield Sample("ggH", "SM", production)
  elif productionmode == "WH":
    yield Sample("WplusH", "SM", production, "POWHEG")
    yield Sample("WminusH", "SM", production, "POWHEG")
  elif productionmode == "ttH":
    yield Sample(productionmode, "SM", "Hff0+", production, "POWHEG")
  else:
    yield Sample(productionmode, "SM", production, "POWHEG")

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
      for color, sample in enumerate(samples(productionmode), start=1):
        if color == 5: color = ROOT.kYellow+3
        hname = "{}{}{}{}{}{}".format(discname, sample.productionmode, sample.hypothesis, sample.hffhypothesis, sample.alternategenerator, analysis)
        wt = "({}>-998)*({})*({}=={})".format(discname, sample.weightname(), categoryname, categoryid)
        tfiles[sample.withdiscriminantsfile()].candTree.Draw("{}>>{}({},{},{})".format(discname, hname, bins, min, max), wt, "hist")
        h = cache[analysis,discname,sample] = getattr(ROOT, hname)
        print h
        h.Scale(1/h.Integral())
        h.SetLineColor(color)
        hstack.Add(h)
        legend.AddEntry(h, str(sample).replace(" "+str(sample.production), ""), "l")

      c = ROOT.TCanvas()
      hstack.Draw("hist nostack")
      hstack.GetXaxis().SetTitle(title)
      legend.Draw()

      directory = os.path.join(config.plotsbasedir, "templateprojections", "compare_POWHEG_JHUGen", str(analysis), str(productionmode))
      try:
        os.makedirs(directory)
      except OSError:
        pass
      for ext in "png eps root pdf".split():
        c.SaveAs(os.path.join(directory, "{}_{}.{}".format(discname, productionmode, ext)))
