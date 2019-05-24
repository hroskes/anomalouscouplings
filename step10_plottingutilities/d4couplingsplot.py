#!/usr/bin/env python

from helperstuff.templates import Template
from helperstuff.enums import *
from helperstuff.utilities import PlotCopier, mkdir_p
import ROOT
import style

class HistHolder(object):
  def __init__(self, h=None):
    self.__h = h
  def __iadd__(self, other):
    if self.__h is None: self.__h = other
    else: self.__h.Add(other)
    return self
  @property
  def hist(self):
    return self.__h

def getproj(category, hypothesis):
  h3D = HistHolder()
  normalization = HistHolder()

  productionmodes, productionmodes_normalization = {
    Category("Untagged"): ({ProductionMode("ggH")}, {ProductionMode("ttH"), ProductionMode("bbH")}),
    Category("VBFtagged"): ({ProductionMode("VBF")}, {ProductionMode("ZH"), ProductionMode("WH")}),
    Category("VHHadrtagged"): ({ProductionMode("ZH"), ProductionMode("WH")}, {ProductionMode("VBF")}),
  }[category]

  if hypothesis is None:
    productionmodes_normalization = set()
    productionmodes = {ProductionMode("qqZZ")}

  for productionmode in productionmodes|productionmodes_normalization:
    if hypothesis == "L1Zg" and productionmode == "WH": continue
    for channel in channels:
      for category2 in categories:
        t = Template(productionmode, hypothesis, "GEN_181119", category, channel, "fa3fa2fL1fL1Zg").gettemplate()
        normalization += t
        if category2 == category and productionmode in productionmodes: h3D += t

  h3D = h3D.hist
  normalization = normalization.hist

  h3D.Scale(1/normalization.Integral())

  return h3D.ProjectionX()

def plot(category, plotcopier=ROOT):
  category = Category(category)

  SM = getproj(category, "0+")
  a3 = getproj(category, "0-")
  a2 = getproj(category, "0h+")
  L1 = getproj(category, "L1")
  L1Zg = getproj(category, "L1Zg")
#  bkg = getproj(category, None)

  hstack = ROOT.THStack("hstack", "hstack")
  hstack.Add(SM)
  hstack.Add(a3)
  hstack.Add(a2)
  hstack.Add(L1)
  hstack.Add(L1Zg)
#  hstack.Add(bkg)

  SM.SetLineColor(2)
  a3.SetLineColor(4)
  a2.SetLineColor(ROOT.kGreen+3)
  L1.SetLineColor(6)
  L1Zg.SetLineColor(ROOT.kViolet-1)
#  bkg.SetLineColor(1)

  legend = ROOT.TLegend(.6, .5, .9, .9)
  legend.SetBorderSize(0)
  legend.SetFillStyle(0)
  legend.AddEntry(SM, "SM", "l")
  legend.AddEntry(a3, "0^{-}", "l")
  legend.AddEntry(a2, "0_{h}^{+}", "l")
  legend.AddEntry(L1, "#Lambda_{1}", "l")
  legend.AddEntry(L1Zg, "#Lambda_{1}^{Z#gamma}", "l")
#  legend.AddEntry(bkg, "qq#rightarrowZZ", "l")

  c = plotcopier.TCanvas()

  hstack.Draw("hist nostack")
  legend.Draw()

  folder = os.path.join(config.plotsbasedir, "templateprojections", "D_4couplings")
  mkdir_p(folder)

  for ext in "png", "root", "pdf":
    c.SaveAs(os.path.join(folder, str(category)+"."+ext))

with PlotCopier() as pc:
  plot("Untagged", pc)
  plot("VBFtagged", pc)
  plot("VHHadrtagged", pc)
