#!/usr/bin/env python

from helperstuff.discriminants import discriminant
from helperstuff.filemanager import tfiles
from helperstuff.enums import *
from helperstuff.samples import *
from helperstuff import style
import ROOT

discriminants = [discriminant(_) for _ in "D_0minus_VBF", "D_CP_VBF", "D_0minus_ZH_hadronic", "D_CP_ZH_hadronic"]

assert len(config.productionsforcombine) == 1
production = config.productionsforcombine[0]

def samples(productionmode):
  for hypothesis in hffhypotheses:
    yield Sample(productionmode, hypothesis, production)
  if productionmode == "HJJ":
    yield Sample("ggH", "SM", production)

cache = {}

for discname, title, bins, min, max in discriminants:
  for productionmode in "HJJ", "ttH":
    hstack = ROOT.THStack(productionmode, title)
    legend = ROOT.TLegend(.6, .7, .9, .9)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    for color, sample in enumerate(samples(productionmode), start=1):
      hname = "{}{}{}".format(discname, sample.productionmode, sample.hypothesis)
      wt = "({}>-998)*{}".format(discname, sample.weightname())
      tfiles[sample.withdiscriminantsfile()].candTree.Draw("{}>>{}({},{},{})".format(discname, hname, bins, min, max), wt, "hist")
      h = cache[discname,sample] = getattr(ROOT, hname)
      h.Scale(1/h.Integral())
      h.SetLineColor(color)
      hstack.Add(h)
      legend.AddEntry(h, str(sample.reweightingsample), "l")

    c = ROOT.TCanvas()
    hstack.Draw("hist nostack")
    hstack.GetXaxis().SetTitle(title)
    legend.Draw()

    for ext in "png eps root pdf".split():
      c.SaveAs(os.path.join(config.plotsbasedir, "TEST", "{}_{}.{}".format(discname, productionmode, ext)))
