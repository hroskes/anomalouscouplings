#!/usr/bin/env python

from helperstuff import style
from helperstuff.discriminants import discriminant
from helperstuff.enums import *
from helperstuff.samples import *
from helperstuff.utilities import tfiles
import ROOT

discriminants = [discriminant(_) for _ in "D_0minus_VBF", "D_CP_VBF", "D_0minus_HadZH", "D_CP_HadZH"]

assert len(config.productionsforcombine) == 1
production = config.productionsforcombine[0]

def category(discname):
    if "VBF" in discname: result = Category("VBFtagged").idnumbers
    if "ZH" in discname: result = Category("VHHadrtagged").idnumbers
    assert len(result) == 1
    return result[0]

def samples(productionmode):
  for hypothesis in hffhypotheses:
    yield Sample(productionmode, hypothesis, production)
  if productionmode == "HJJ":
    yield Sample("ggH", "SM", production)
  elif productionmode == "WH":
    yield Sample("WplusH", "SM", production, "POWHEG")
    yield Sample("WminusH", "SM", production, "POWHEG")
  else:
    yield Sample(productionmode, "SM", production, "POWHEG")

cache = {}

for discname, title, bins, min, max in discriminants:
  for productionmode in "HJJ", "ttH", "VBF", "ZH", "WH":
    hstack = ROOT.THStack(productionmode, title)
    legend = ROOT.TLegend(.6, .7, .9, .9)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    for color, sample in enumerate(samples(productionmode), start=1):
      hname = "{}{}{}{}".format(discname, sample.productionmode, sample.hypothesis, sample.alternategenerator)
      wt = "({}>-998)*({})*(category=={})".format(discname, sample.weightname(), category(discname))
      tfiles[sample.withdiscriminantsfile()].candTree.Draw("{}>>{}({},{},{})".format(discname, hname, bins, min, max), wt, "hist")
      h = cache[discname,sample] = getattr(ROOT, hname)
      print h
      h.Scale(1/h.Integral())
      h.SetLineColor(color)
      hstack.Add(h)
      legend.AddEntry(h, str(sample).replace(" "+str(sample.production), ""), "l")

    c = ROOT.TCanvas()
    hstack.Draw("hist nostack")
    hstack.GetXaxis().SetTitle(title)
    legend.Draw()

    directory = os.path.join(config.plotsbasedir, "templateprojections", "compare_POWHEG_JHUGen", str(productionmode))
    try:
      os.makedirs(directory)
    except OSError:
      pass
    for ext in "png eps root pdf".split():
      c.SaveAs(os.path.join(directory, "{}_{}.{}".format(discname, productionmode, ext)))
