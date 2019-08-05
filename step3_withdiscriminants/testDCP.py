#!/usr/bin/env python

import os

from collections import namedtuple

import ROOT

from helperstuff import config, style
from helperstuff.plotfromtree import plotfromtree
from helperstuff.samples import ReweightingSamplePlus
from helperstuff.utilities import PlotCopier

commonkwargs = dict(
  disc="D_CP_VBF_2bins",
  cut="D_bkg_VBFdecay>0.2",
  category="VBFtagged",
  normalizeto1=True,
  analysis="fa3fa2fL1fL1Zg_morecategories",
  production="190703_2017",
  masscut=True
)

neededkwargs = [
  "reweightfrom",
  "color",
  "hname",
]

Line = namedtuple("Line", neededkwargs)

lines = [
  Line(ReweightingSamplePlus("VBF", "0+"), 4, "VBFSM"),
  Line(ReweightingSamplePlus("VBF", "fa3prod0.5"), ROOT.kAzure-4, "VBFmix"),
  Line(ReweightingSamplePlus("HJJ", "0+", "Hff0+", "ext2"), ROOT.kGreen+3, "ggHSM"),
  Line(ReweightingSamplePlus("HJJ", "0+", "fCP0.5", "ext2"), ROOT.kGreen-3, "ggHmix"),
  Line(ReweightingSamplePlus("ttH", "0+", "Hff0+"), ROOT.kMagenta+3, "ttHSM"),
  Line(ReweightingSamplePlus("ttH", "0+", "fCP0.5"), ROOT.kMagenta-4, "ttHmix"),
]

cache = []

hstack = ROOT.THStack("hstack", "")

for line in lines:
  kwargs = commonkwargs.copy()
  kwargs.update(**line._asdict())
  h = plotfromtree(**kwargs)
  cache.append(h)
  hstack.Add(h)

l = ROOT.TLegend(.6, .5, .9, .9)
l.SetBorderSize(0)
l.SetFillStyle(0)

with PlotCopier() as pc:
  c = pc.TCanvas()
  hstack.Draw("hist nostack")
  l.Draw()
  for ext in "png root pdf C".split():
    c.SaveAs(os.path.join(config.plotsbasedir, "TEST", "DCPVBF."+ext))
