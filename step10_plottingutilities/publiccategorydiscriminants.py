#!/usr/bin/env python

import argparse
import glob
import itertools
import os
import re

if __name__ == "__main__":
  p = argparse.ArgumentParser()
  p.add_argument("--analysis", choices="fa3 fa2 fL1 fL1Zg".split())
  p.add_argument("--VBForVH", choices="VBF VH".split())
  p.add_argument("--cut", choices="decay proddec bycategory".split())
  args = p.parse_args()

import ROOT

from helperstuff import config, stylefunctions as style
from helperstuff.combinehelpers import Luminosity
from helperstuff.enums import Analysis, flavors
from helperstuff.samples import ReweightingSample, Sample
from helperstuff.utilities import PlotCopier, TFile

class HistAndNormalization(object):
  def __init__(self, h, hnormalization):
    self.__h = h
    self.__hnormalization = hnormalization
  def Scale(self, *args):
    self.__hnormalization.Scale(*args)
    return self.__h.Scale(*args)
  def Integral(self, *args):
    return self.__hnormalization.Integral(*args)
  def Add(self, other):
    self.__h.Add(other.__h)
    self.__hnormalization.Add(other.__hnormalization)
  @property
  def h(self):
    return self.__h
  def __getattr__(self, attr):
    if re.match("[SG]et(Line|Marker|Fill)(Size|Color|Style|Width)", attr): return getattr(self.__h, attr)
    raise AttributeError("'{}' object has no attribute '{}'".format(type(self).__name__, attr))

def histogram(analysis, VBForVH, *samples, **kwargs):
  reweightto = kwargs.pop("reweightto", None)
  whichdbkg = kwargs.pop("whichdbkg", "decay")
  assert not kwargs, kwargs

  productionmode = {s.productionmode for s in samples}
  assert len(productionmode) == 1, productionmode
  productionmode = productionmode.pop()

  filenames = [sample.withdiscriminantsfile() for sample in samples]

  histname = "".join(os.path.basename(filename).replace(".root", "") for filename in filenames) + ("reweighted" if reweightto else "")
  h = ROOT.TH1F(histname, "", 10, 0, 1)
  hnormalization = ROOT.TH1F(histname+"_normalization", "", 1, 0, 1)

  analysis = Analysis(analysis)
  discriminantnames = {
    "VBF": ("D_2jet_0plus",),
    "VH": ("D_HadZH_0plus", "D_HadWH_0plus"),
  }[VBForVH]
  discriminantnames += tuple(
    _.replace(
      "0plus", {
        Analysis("fa3"): "0minus",
        Analysis("fa2"): "a2",
        Analysis("fL1"): "L1",
        Analysis("fL1Zg"): "L1Zg",
      }[fai]
    ) for fai in analysis.fais for _ in discriminantnames
  )

  categorizationname = {
    Analysis("fa3"): "category_0P_or_0M",
    Analysis("fa2"): "category_0P_or_a2",
    Analysis("fL1"): "category_0P_or_L1",
    Analysis("fL1Zg"): "category_0P_or_L1Zg",
  }[analysis]
  dbkgname = {
    "VBF": "D_bkg_VBFdecay",
    "VH": "D_bkg_HadVHdecay",
  }[VBForVH]

  xformula = reduce(lambda x, y: "max("+x+", "+y+")", discriminantnames)
  mastercutformula = "105 < ZZMass && ZZMass < 140"
  cutformula = {
    "VBF": "nExtraLep==0 && (((nCleanedJetsPt30==2||nCleanedJetsPt30==3)&&nCleanedJetsPt30BTagged_bTagSF<=1)||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0))",
    "VH": "nExtraLep==0 && (nCleanedJetsPt30==2||nCleanedJetsPt30==3||(nCleanedJetsPt30>=4&&nCleanedJetsPt30BTagged_bTagSF==0))"
  }[VBForVH] + {
    "bycategory": " && (({cat} == 2 && D_bkg_VBFdecay > 0.5) || ({cat} == 4 && D_bkg_HadVHdecay > 0.5) || ({cat} != 2 && {cat} != 4 && D_bkg > 0.5))",
    "proddec": " && {dbkgname} > 0.5",
    "decay": " && D_bkg > 0.5",
  }[whichdbkg].format(cat=categorizationname, dbkgname=dbkgname)
  if productionmode == "data":
    weightformula = "1"
  elif productionmode == "ZX":
    weightformula = "MC_weight_ZX"
  elif reweightto is not None:
    weightformula = reweightto.weightname()
  else:
    weightformula = "overallEventWeight"

  for filename in filenames:
    with TFile(filename) as f:
      t = f.candTree

      xtreeformula = ROOT.TTreeFormula("x", xformula, t)
      mastercuttreeformula = ROOT.TTreeFormula("mastercut", mastercutformula, t)
      cuttreeformula = ROOT.TTreeFormula("cut", cutformula, t)
      weighttreeformula = ROOT.TTreeFormula("weight", weightformula, t)

      t.SetBranchStatus("*", 0)
      for formula in xformula, mastercutformula, cutformula, weightformula:
        for branch in re.findall(r"\b[a-zA-Z_][a-zA-Z0-9_]+\b", formula):
          if branch == "max": continue
          t.SetBranchStatus(branch, 1)
          try:
            t.GetEntry(0)
            getattr(t, branch)
          except AttributeError:
            raise ValueError("Bad branch "+branch+" in "+str(productionmode))

      size = t.GetEntries()
      print filename
      for i, entry in enumerate(t, start=1):
        if mastercuttreeformula.EvalInstance():
          if cuttreeformula.EvalInstance():
            h.Fill(xtreeformula.EvalInstance(), weighttreeformula.EvalInstance())
          hnormalization.Fill(0, weighttreeformula.EvalInstance())
        if i % 10000 == 0 or i == size:
          print i, "/", size
          #break

  return HistAndNormalization(h, hnormalization)

def publiccategorydiscriminants(analysis, VBForVH, plotcopier=ROOT, whichdbkg="decay"):
  analysis = Analysis(analysis)

  hstack = ROOT.THStack()

  data = histogram(analysis, VBForVH, Sample("data", "180721"), Sample("data", "180722"), whichdbkg=whichdbkg)

  ZX = histogram(analysis, VBForVH, Sample("ZX", "180721"), Sample("ZX", "180722"), whichdbkg=whichdbkg)
  qqZZ = histogram(analysis, VBForVH, Sample("qqZZ", "180721"), Sample("qqZZ", "180722"), whichdbkg=whichdbkg)
  ggZZ = histogram(analysis, VBForVH, *(Sample("ggZZ", year, flavor) for year in ("180721", "180722") for flavor in flavors), whichdbkg=whichdbkg)
  VBFbkg = histogram(analysis, VBForVH, *(Sample("VBFbkg", year, flavor) for year in ("180721",) for flavor in ("2e2mu", "4e", "4mu")), whichdbkg=whichdbkg)

  SMhypothesis = "0+"
  BSMhypothesis = {
    Analysis("fa3"): "0-",
    Analysis("fa2"): "0h+",
    Analysis("fL1"): "L1",
    Analysis("fL1Zg"): "L1Zg",
  }[analysis]

  VBFSM  = histogram(analysis, VBForVH, Sample("VBF", SMhypothesis,  "180721"), Sample("VBF", SMhypothesis,  "180722"), whichdbkg=whichdbkg)
  ZHSM   = histogram(analysis, VBForVH, Sample("ZH",  SMhypothesis,  "180721"), Sample("ZH",  SMhypothesis,  "180722"), whichdbkg=whichdbkg)
  WHSM   = histogram(analysis, VBForVH, Sample("WH",  SMhypothesis,  "180721"), Sample("WH",  SMhypothesis,  "180722"), whichdbkg=whichdbkg)
  ggHSM  = histogram(analysis, VBForVH, Sample("ggH", SMhypothesis,  "180721"), Sample("ggH", SMhypothesis,  "180722"), whichdbkg=whichdbkg)
  ttHSM  = histogram(analysis, VBForVH, Sample("ttH", "Hff0+", SMhypothesis,  "180721"), Sample("ttH", "Hff0+", SMhypothesis,  "180722"), whichdbkg=whichdbkg)
  bbHSM  = histogram(analysis, VBForVH, Sample("bbH", SMhypothesis,  "180721"), Sample("bbH", SMhypothesis,  "180722"), whichdbkg=whichdbkg)

  if BSMhypothesis == "L1Zg":
    VBFBSM = histogram(analysis, VBForVH, Sample("VBF", "L1", "180721"), Sample("VBF", "L1", "180722"), reweightto=ReweightingSample("VBF", BSMhypothesis), whichdbkg=whichdbkg)
    ZHBSM  = histogram(analysis, VBForVH, Sample("ZH",  "L1", "180721"), Sample("ZH",  "L1", "180722"), reweightto=ReweightingSample("ZH", BSMhypothesis), whichdbkg=whichdbkg)
    WHBSM  = histogram(analysis, VBForVH, Sample("WH",  "L1", "180721"), Sample("WH",  "L1", "180722"), reweightto=ReweightingSample("WH", BSMhypothesis), whichdbkg=whichdbkg)
    ggHBSM = histogram(analysis, VBForVH, Sample("ggH", "0+", "180721"), Sample("ggH", "0+", "180722"), reweightto=ReweightingSample("ggH", BSMhypothesis), whichdbkg=whichdbkg)
  else:
    VBFBSM = histogram(analysis, VBForVH, Sample("VBF", BSMhypothesis, "180721"), Sample("VBF", BSMhypothesis, "180722"), whichdbkg=whichdbkg)
    ZHBSM  = histogram(analysis, VBForVH, Sample("ZH",  BSMhypothesis, "180721"), Sample("ZH",  BSMhypothesis, "180722"), whichdbkg=whichdbkg)
    WHBSM  = histogram(analysis, VBForVH, Sample("WH",  BSMhypothesis, "180721"), Sample("WH",  BSMhypothesis, "180722"), whichdbkg=whichdbkg)
    ggHBSM = histogram(analysis, VBForVH, Sample("ggH", BSMhypothesis, "180721"), Sample("ggH", BSMhypothesis, "180722"), whichdbkg=whichdbkg)
  ttHBSM = histogram(analysis, VBForVH, Sample("ttH", "Hff0+", "0+", "180721"), Sample("ttH", "Hff0+", "0+", "180722"), reweightto=ReweightingSample("ttH", "Hff0+", BSMhypothesis), whichdbkg=whichdbkg)
  bbHBSM = histogram(analysis, VBForVH, Sample("bbH", "0+", "180721"), Sample("bbH", "0+", "180722"), reweightto=ReweightingSample("bbH", BSMhypothesis), whichdbkg=whichdbkg)

  VVHSM = VBFSM
  VVHSM.Add(ZHSM)
  VVHSM.Add(WHSM)
  ffHSM = ggHSM
  ffHSM.Add(ttHSM)
  ffHSM.Add(bbHSM)

  VVHBSM = VBFBSM
  VVHBSM.Add(ZHBSM)
  VVHBSM.Add(WHBSM)
  ffHBSM = ggHBSM
  ffHBSM.Add(ttHBSM)
  ffHBSM.Add(bbHBSM)

  ZZ = qqZZ
  ZZ.Add(ggZZ)
  ZZ.Add(VBFbkg)

  ZX.SetLineColor(1)
  ZX.SetLineWidth(2)
  ZX.SetFillColor(ROOT.TColor.GetColor("#669966"))
  ZX.SetFillStyle(1001)

  ZZ.SetLineColor(1)
  ZZ.SetLineWidth(2)
  ZZ.SetFillColor(ROOT.kAzure-9)
  ZZ.SetFillStyle(1001)

  VVHSM.SetLineColor(4)
  VVHSM.SetLineWidth(2)
  VVHSM.SetFillStyle(0)
  ffHSM.SetLineColor(ROOT.kOrange+10)
  ffHSM.SetLineWidth(2)
  ffHSM.SetFillStyle(0)

  VVHBSM.SetLineColor(4)
  VVHBSM.SetLineWidth(2)
  VVHBSM.SetLineStyle(2)
  VVHBSM.SetFillStyle(0)
  ffHBSM.SetLineColor(ROOT.kOrange+10)
  ffHBSM.SetLineWidth(2)
  ffHBSM.SetLineStyle(2)
  ffHBSM.SetFillStyle(0)

  histogramorder = ffHSM, VVHSM, ffHBSM, VVHBSM, ZZ, ZX
  normalization = {h: 0 for h in histogramorder+(data,)}

  nfiles = 0
  g = glob.glob(os.path.join(config.plotsbasedir, "templateprojections/niceplots/fullrange", str(analysis), "*tagged", "D_bkg*.root"))
  for filename in g:
    nfiles += 1
    with TFile(filename) as f:
      for i, h in enumerate(histogramorder):
        proj = f.cprojections.GetListOfPrimitives()[1].GetHists()[i]
        try:
          assert proj.GetMarkerColor() == h.GetMarkerColor(), (h.GetMarkerColor(), proj.GetMarkerColor())
          assert proj.GetMarkerStyle() == h.GetMarkerStyle(), (h.GetMarkerStyle(), proj.GetMarkerStyle())
          assert proj.GetMarkerSize() == h.GetMarkerSize(), (h.GetMarkerSize(), proj.GetMarkerSize())
          assert proj.GetLineColor() == h.GetLineColor(), (h.GetLineColor(), proj.GetLineColor())
          assert proj.GetLineStyle() == h.GetLineStyle(), (h.GetLineStyle(), proj.GetLineStyle())
          assert proj.GetLineWidth() == h.GetLineWidth(), (h.GetLineWidth(), proj.GetLineWidth())
          assert proj.GetFillColor() == h.GetFillColor(), (h.GetFillColor(), proj.GetFillColor())
          assert proj.GetFillStyle() == h.GetFillStyle(), (h.GetFillStyle(), proj.GetFillStyle())
        except AssertionError:
          print i
          f.cprojections.GetListOfPrimitives()[1].ls()
          raise
        normalization[h] += proj.Integral()

      datagraph = f.cprojections.GetListOfPrimitives()[6]
      normalization[data] += sum(y for _, y in itertools.izip(xrange(datagraph.GetN()), datagraph.GetY()))
  assert nfiles == 3, g

  normalization[ffHSM] -= normalization[VVHSM]
  normalization[ffHBSM] -= normalization[VVHBSM]
  normalization[VVHSM] -= normalization[ZZ]
  normalization[VVHBSM] -= normalization[ZZ]
  normalization[ZZ] -= normalization[ZX]

  for h, integral in normalization.iteritems():
    if h == data:
      assert h.Integral() == integral, (h.Integral(), integral)
    else:
      h.Scale(integral / h.Integral())

  data = style.asymmerrorsfromhistogram(data.h, showemptyerrors=False)
  data.SetLineColor(1)
  data.SetMarkerColor(1)
  data.SetLineStyle(1)
  data.SetLineWidth(1)
  data.SetMarkerStyle(20)
  data.SetMarkerSize(1.2)

  ZZ.Add(ZX)
  VVHSM.Add(ZZ)
  VVHBSM.Add(ZZ)

  totalSM = ffHSM
  totalSM.Add(VVHSM)
  totalBSM = ffHBSM
  totalBSM.Add(VVHBSM)

  hstack = ROOT.THStack()
  hstack.Add(ZZ.h, "hist")
  hstack.Add(ZX.h, "hist")
  hstack.Add(VVHSM.h, "hist")
  hstack.Add(VVHBSM.h, "hist")
  hstack.Add(ffHSM.h, "hist")
  hstack.Add(ffHBSM.h, "hist")

  ymax = style.ymax((hstack, "nostack"), (data, "P"))
  if VBForVH == "VBF": ymax = max(ymax, 9.20351838988)
  hstack.SetMaximum(ymax)

  c = plotcopier.TCanvas("c", "", 8, 30, 800, 800)
  style.applycanvasstyle(c)
  hstack.Draw("nostack")

  BSMname = {
    Analysis("fa3"): "0-",
    Analysis("fa2"): "0h+",
    Analysis("fL1"): "#Lambda1",
    Analysis("fL1Zg"): "#Lambda1Z#gamma",
  }[analysis]
  discname = {
    "VBF": "max#left(D_{{2jet}}^{{VBF}}, D_{{2jet}}^{{VBF,{BSMname}}}#right)",
    "VH":  "max#left(D_{{2jet}}^{{WH}}, D_{{2jet}}^{{WH,{BSMname}}}, D_{{2jet}}^{{ZH}}, D_{{2jet}}^{{ZH,{BSMname}}}#right)",
  }[VBForVH].format(BSMname=BSMname).replace(", D_{2jet}^{WH,#Lambda1Z#gamma}", "")
  hstack.GetXaxis().SetTitle(discname)
  hstack.GetYaxis().SetTitle("Events / bin")
  style.applyaxesstyle(hstack)
  hstack.GetXaxis().SetTitleSize(0.05)
  hstack.GetXaxis().SetTitleOffset(1.1)

  hstack.GetXaxis().CenterTitle()
  hstack.GetYaxis().CenterTitle()

  l = ROOT.TLegend(.6, .55, .9, .9)
  style.applylegendstyle(l)
  l.AddEntry(data, "Observed", "lp")
  l.AddEntry(totalSM.h, "Total SM", "l")
  l.AddEntry(VVHSM.h, "VBF+VH SM", "f")
  l.AddEntry(totalBSM.h, "Total {} = 1".format(analysis.title()), "l")
  l.AddEntry(VVHBSM.h, "VBF+VH {} = 1".format(analysis.title()), "f")
  l.AddEntry(ZZ.h, "ZZ/Z#gamma*", "f")
  l.AddEntry(ZX.h, "Z+X", "f")

  line = ROOT.TLine(0.5, 0, 0.5, ymax)
  line.SetLineWidth(4)
  line.SetLineColor(ROOT.kViolet-1)
  line.SetLineStyle(9)
  line.Draw()

  data.Draw("PE")

  style.CMS("" if analysis == "fa3" else "Supplementary", sum(float(Luminosity("fordata", production)) for production in config.productionsforcombine))

  l.Draw()

  for ext in "png root pdf eps C".split():
    c.SaveAs(os.path.join(config.plotsbasedir, "templateprojections/niceplots/fullrange", str(analysis), "D_2jet_"+VBForVH+"cut"+whichdbkg+"."+ext))

if __name__ == "__main__":
  with PlotCopier() as pc:
    for analysis in "fa3", "fa2", "fL1", "fL1Zg":
      for VBForVH in "VBF", "VH":
        for cut in "decay", "proddec", "bycategory":
          if analysis != args.analysis is not None: continue
          if VBForVH != args.VBForVH is not None: continue
          if cut != args.cut is not None: continue
          publiccategorydiscriminants(analysis, VBForVH, pc, whichdbkg=cut)
