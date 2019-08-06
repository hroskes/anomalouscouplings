#!/usr/bin/env python

import os

from collections import namedtuple

import ROOT

from helperstuff import config, style
from helperstuff.plotfromtree import plotfromtree
from helperstuff.samples import ReweightingSamplePlus
from helperstuff.utilities import PlotCopier

CommonKwargs = namedtuple("CommonKwargs", [
  "disc",
  "cut",
  "category",
  "normalizeto1",
  "analysis",
  "production",
  "masscut",
])

neededkwargs = [
  "reweightto",
  "color",
  "hname",
  "reweightfrom",
]

class Line(namedtuple("Line", neededkwargs)):
  def __new__(cls, reweightto, color, hname, reweightfrom=None):
    if reweightfrom is None: reweightfrom=reweightto
    return super(Line, cls).__new__(cls, reweightto, color, hname, reweightfrom)

def run(commonkwargs, lines, saveas, pc=ROOT):

  cache = []

  hstack = ROOT.THStack("hstack", "")

  l = ROOT.TLegend(.6, .5, .9, .9)
  l.SetBorderSize(0)
  l.SetFillStyle(0)

  for line in lines:
    kwargs = commonkwargs._asdict()
    kwargs.update(**line._asdict())
    kwargs["hname"] += kwargs["disc"]
    h = plotfromtree(**kwargs)
    cache.append(h)
    hstack.Add(h)
    l.AddEntry(h, line.hname, "l")

  c = pc.TCanvas()
  hstack.Draw("hist nostack")
  l.Draw()
  for ext in "png root pdf C".split():
    c.SaveAs(os.path.join(config.plotsbasedir, "TEST", saveas+"."+ext))

with PlotCopier() as pc:
  DCPVBFkwargs = CommonKwargs(
    disc="D_CP_VBF_2bins",
    cut="D_bkg_VBFdecay>0.2",
    category="VBFtagged",
    normalizeto1=True,
    analysis="fa3fa2fL1fL1Zg_morecategories",
    production="190703_2017",
    masscut=True,
  )

  D0mVBFkwargs = CommonKwargs(
    disc="D_0minus_VBFdecay_3bins",
    cut="D_bkg_VBFdecay>0.2",
    category="VBFtagged",
    normalizeto1=True,
    analysis="fa3fa2fL1fL1Zg_morecategories",
    production="190703_2017",
    masscut=True,
  )

  DCPVHkwargs = CommonKwargs(
    disc="D_CP_HadVH_2bins",
    cut="D_bkg_HadVHdecay>0.2",
    category="VHHadrtagged",
    normalizeto1=True,
    analysis="fa3fa2fL1fL1Zg_morecategories",
    production="190703_2017",
    masscut=True,
  )

  D0mVHkwargs = CommonKwargs(
    disc="D_0minus_HadVHdecay_3bins",
    cut="D_bkg_HadVHdecay>0.2",
    category="VHHadrtagged",
    normalizeto1=True,
    analysis="fa3fa2fL1fL1Zg_morecategories",
    production="190703_2017",
    masscut=True,
  )

  VBFmixlines = [
    Line(ReweightingSamplePlus("VBF", "0+"), 4, "VBFSM"),
    Line(ReweightingSamplePlus("VBF", "fa3prod0.5"), ROOT.kAzure-4, "VBFmix"),
    Line(ReweightingSamplePlus("ggH", "0+", "Hff0+", "MC@NLO"), ROOT.kGreen+3, "ggHSM"),
    Line(ReweightingSamplePlus("ggH", "0+", "fCP0.5", "MC@NLO"), ROOT.kGreen-3, "ggHmix"),
    Line(ReweightingSamplePlus("ttH", "0+", "Hff0+"), ROOT.kMagenta+3, "ttHSM"),
    Line(ReweightingSamplePlus("ttH", "0+", "fCP0.5"), ROOT.kMagenta-4, "ttHmix"),
  ]

  VBFpurelines = [
    Line(ReweightingSamplePlus("VBF", "0+"), 4, "VBFSM"),
    Line(ReweightingSamplePlus("VBF", "0-"), ROOT.kAzure-4, "VBFPS", reweightfrom=ReweightingSamplePlus("VBF", "fa3prod0.5")),
    Line(ReweightingSamplePlus("ggH", "0+", "Hff0+", "MC@NLO"), ROOT.kGreen+3, "ggHSM"),
    Line(ReweightingSamplePlus("ggH", "0+", "Hff0-", "MC@NLO"), ROOT.kGreen-3, "ggHPS"),
    Line(ReweightingSamplePlus("ttH", "0+", "Hff0+"), ROOT.kMagenta+3, "ttHSM"),
    Line(ReweightingSamplePlus("ttH", "0+", "Hff0-", "ext1"), ROOT.kMagenta-4, "ttHPS"),
  ]

  VHmixlines = [
    Line(ReweightingSamplePlus("ZH", "0+", "ext1"), 4, "ZHSM"),
    Line(ReweightingSamplePlus("ZH", "fa3prod0.5"), ROOT.kAzure-4, "ZHmix"),
    Line(ReweightingSamplePlus("ggH", "0+", "Hff0+", "MC@NLO"), ROOT.kGreen+3, "ggHSM"),
    Line(ReweightingSamplePlus("ggH", "0+", "fCP0.5", "MC@NLO"), ROOT.kGreen-3, "ggHmix"),
    Line(ReweightingSamplePlus("ttH", "0+", "Hff0+"), ROOT.kMagenta+3, "ttHSM"),
    Line(ReweightingSamplePlus("ttH", "0+", "fCP0.5"), ROOT.kMagenta-4, "ttHmix"),
  ]

  VHpurelines = [
    Line(ReweightingSamplePlus("ZH", "0+", "ext1"), 4, "ZHSM"),
    Line(ReweightingSamplePlus("ZH", "0-"), ROOT.kAzure-4, "ZHPS"),
    Line(ReweightingSamplePlus("ggH", "0+", "Hff0+", "MC@NLO"), ROOT.kGreen+3, "ggHSM"),
    Line(ReweightingSamplePlus("ggH", "0+", "Hff0-", "MC@NLO"), ROOT.kGreen-3, "ggHPS"),
    Line(ReweightingSamplePlus("ttH", "0+", "Hff0+"), ROOT.kMagenta+3, "ttHSM"),
    Line(ReweightingSamplePlus("ttH", "0+", "Hff0-", "ext1"), ROOT.kMagenta-4, "ttHPS"),
  ]

  run(DCPVBFkwargs, VBFmixlines, "D_CP_VBF", pc)
  run(D0mVBFkwargs, VBFpurelines, "D_0minus_VBFdecay", pc)
  run(DCPVHkwargs, VHmixlines, "D_CP_VH", pc)
  run(D0mVHkwargs, VHpurelines, "D_0minus_VHdecay", pc)
