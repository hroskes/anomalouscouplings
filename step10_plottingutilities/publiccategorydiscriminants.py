#!/usr/bin/env python

import argparse
import glob
import itertools
import os
import re

if __name__ == "__main__":
  p = argparse.ArgumentParser()
  p.add_argument("--submit", action="store_true")
  p.add_argument("--analysis", choices="fa3 fa2 fL1 fL1Zg".split())
  p.add_argument("--VBForVH", choices="VBF VH".split())
  p.add_argument("--cut", choices="decay proddec bycategory m4l".split())
  args = p.parse_args()

  if args.submit:
    for analysis in "fa3 fa2 fL1 fL1Zg".split():
      for VBForVH in "VBF VH".split():
        for cut in "decay proddec bycategory m4l".split():
          os.system("""echo "cd $CMSSW_BASE && eval "'$(scram ru -sh)'" && cd $(pwd) && ./publiccategorydiscriminants.py --a {}  --V {} --c {}" | bsub -q cmscaf1nd""".format(analysis, VBForVH, cut))
    import sys
    sys.exit(0)

import ROOT

from helperstuff import config, stylefunctions as style
from helperstuff.combinehelpers import getrate, Luminosity
from helperstuff.enums import Analysis, categories, channels, flavors, Production
from helperstuff.samples import ReweightingSample, ReweightingSamplePlus, ReweightingSampleWithFlavor, Sample
from helperstuff.utilities import PlotCopier, TFile

productions = {Production(p) for p in ("180721", "180722")}

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
  def Clone(self, name):
    return HistAndNormalization(self.__h.Clone(name), self.__hnormalization.Clone(name+"_normalization"))
  @property
  def h(self):
    return self.__h
  def __getattr__(self, attr):
    if re.match("[SG]et(Line|Marker|Fill)(Size|Color|Style|Width)", attr): return getattr(self.__h, attr)
    raise AttributeError("'{}' object has no attribute '{}'".format(type(self).__name__, attr))

__scalefactors = {}

def histogram(analysis, VBForVH, *reweightingsamples, **kwargs):
  reweightto = kwargs.pop("reweightto", None)
  whichdbkg = kwargs.pop("whichdbkg", "decay")
  assert not kwargs, kwargs

  productionmode = {s.productionmode for s in reweightingsamples}
  assert len(productionmode) == 1, productionmode
  productionmode = productionmode.pop()

  hypothesis = {s.hypothesis for s in reweightingsamples}
  assert len(hypothesis) == 1, hypothesis
  hypothesis = hypothesis.pop()

  histname = str(reweightingsamples[0]) + ("reweighted" if reweightto else "")
  finalh = ROOT.TH1F(histname, "", 10, 0, 1)
  finalhnormalization = ROOT.TH1F(histname+"_normalization", "", 1, 0, 1)
  finalresult = HistAndNormalization(finalh, finalhnormalization)

  for production in productions:
    if production == "180722" and productionmode == "VBF bkg": continue

    filenames = [Sample(sample, production).withdiscriminantsfile() for sample in reweightingsamples]

    h = ROOT.TH1F("tmp"+str(production), "", 10, 0, 1)
    hnormalization = ROOT.TH1F("tmp_normalization"+str(production), "", 1, 0, 1)
    handnormalization = HistAndNormalization(h, hnormalization)

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
      "m4l": " && 118 < ZZMass && ZZMass < 130"
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

    if hypothesis is None or hypothesis == "SM":
      if productionmode == "data":
        __scalefactors[productionmode, production] = 1
      else:
        __scalefactors[productionmode, production] = sum(getrate(ch, ca, productionmode, "fordata", production, analysis) for ch in channels for ca in categories) / handnormalization.Integral()

    handnormalization.Scale(__scalefactors[productionmode, production])

    finalresult.Add(handnormalization)

  return finalresult

def publiccategorydiscriminants(analysis, VBForVH, plotcopier=ROOT, whichdbkg="decay"):
  analysis = Analysis(analysis)

  hstack = ROOT.THStack()

  data = histogram(analysis, VBForVH, ReweightingSample("data"), whichdbkg=whichdbkg)

  ZX = histogram(analysis, VBForVH, ReweightingSample("ZX"), whichdbkg=whichdbkg)
  qqZZ = histogram(analysis, VBForVH, ReweightingSample("qqZZ"), ReweightingSamplePlus("qqZZ", "ext"), whichdbkg=whichdbkg)
  ggZZ = histogram(analysis, VBForVH, *(ReweightingSampleWithFlavor("ggZZ", flavor) for flavor in flavors), whichdbkg=whichdbkg)
  VBFbkg = histogram(analysis, VBForVH, *(ReweightingSampleWithFlavor("VBFbkg", flavor) for flavor in ("2e2mu", "4e", "4mu")), whichdbkg=whichdbkg)

  SMhypothesis = "0+"
  BSMhypothesis = {
    Analysis("fa3"): "0-",
    Analysis("fa2"): "0h+",
    Analysis("fL1"): "L1",
    Analysis("fL1Zg"): "L1Zg",
  }[analysis]

  VBFSM  = histogram(analysis, VBForVH, ReweightingSample("VBF", SMhypothesis), whichdbkg=whichdbkg)
  ZHSM   = histogram(analysis, VBForVH, ReweightingSample("ZH",  SMhypothesis), whichdbkg=whichdbkg)
  WHSM   = histogram(analysis, VBForVH, ReweightingSample("WH",  SMhypothesis), whichdbkg=whichdbkg)
  ggHSM  = histogram(analysis, VBForVH, ReweightingSample("ggH", SMhypothesis), whichdbkg=whichdbkg)
  ttHSM  = histogram(analysis, VBForVH, ReweightingSample("ttH", "Hff0+", SMhypothesis), whichdbkg=whichdbkg)
  bbHSM  = histogram(analysis, VBForVH, ReweightingSample("bbH", SMhypothesis), whichdbkg=whichdbkg)

  if BSMhypothesis == "L1Zg":
    VBFBSM = histogram(analysis, VBForVH, ReweightingSample("VBF", "L1"), reweightto=ReweightingSample("VBF", BSMhypothesis), whichdbkg=whichdbkg)
    ZHBSM  = histogram(analysis, VBForVH, ReweightingSample("ZH",  "L1"), reweightto=ReweightingSample("ZH", BSMhypothesis), whichdbkg=whichdbkg)
    WHBSM  = histogram(analysis, VBForVH, ReweightingSample("WH",  "L1"), reweightto=ReweightingSample("WH", BSMhypothesis), whichdbkg=whichdbkg)
    ggHBSM = histogram(analysis, VBForVH, ReweightingSample("ggH", "0+"), reweightto=ReweightingSample("ggH", BSMhypothesis), whichdbkg=whichdbkg)
  else:
    VBFBSM = histogram(analysis, VBForVH, ReweightingSample("VBF", BSMhypothesis), whichdbkg=whichdbkg)
    ZHBSM  = histogram(analysis, VBForVH, ReweightingSample("ZH",  BSMhypothesis), whichdbkg=whichdbkg)
    WHBSM  = histogram(analysis, VBForVH, ReweightingSample("WH",  BSMhypothesis), whichdbkg=whichdbkg)
    ggHBSM = histogram(analysis, VBForVH, ReweightingSample("ggH", BSMhypothesis), whichdbkg=whichdbkg)
  ttHBSM = histogram(analysis, VBForVH, ReweightingSample("ttH", "Hff0+", "0+"), reweightto=ReweightingSample("ttH", "Hff0+", BSMhypothesis), whichdbkg=whichdbkg)
  bbHBSM = histogram(analysis, VBForVH, ReweightingSample("bbH", "0+"), reweightto=ReweightingSample("bbH", BSMhypothesis), whichdbkg=whichdbkg)

  scaleVVHBSM =  (VBFSM.Integral() + ZHSM.Integral() + WHSM.Integral()) / (VBFBSM.Integral() + ZHBSM.Integral() + WHBSM.Integral())
  VBFBSM.Scale(scaleVVHBSM)
  ZHBSM.Scale(scaleVVHBSM)
  WHBSM.Scale(scaleVVHBSM)
  scaleffHBSM =  (ggHSM.Integral() + ttHSM.Integral() + bbHSM.Integral()) / (ggHBSM.Integral() + ttHBSM.Integral() + bbHBSM.Integral())
  ggHBSM.Scale(scaleffHBSM)
  ttHBSM.Scale(scaleffHBSM)
  bbHBSM.Scale(scaleffHBSM)

  othersignalSM = ggHSM.Clone("othersignalSM")
  othersignalSM.Add(ttHSM)
  othersignalSM.Add(bbHSM)
  if VBForVH == "VBF":
    targetSM = VBFSM.Clone("targetSM")
    othersignalSM.Add(ZHSM)
    othersignalSM.Add(WHSM)
  elif VBForVH == "VH":
    targetSM = ZHSM.Clone("targetSM")
    targetSM.Add(WHSM)
    othersignalSM.Add(VBFSM)

  othersignalBSM = ggHBSM.Clone("othersignalBSM")
  othersignalBSM.Add(ttHBSM)
  othersignalBSM.Add(bbHBSM)
  if VBForVH == "VBF":
    targetBSM = VBFBSM.Clone("targetBSM")
    othersignalBSM.Add(ZHBSM)
    othersignalBSM.Add(WHBSM)
  elif VBForVH == "VH":
    targetBSM = ZHBSM.Clone("targetBSM")
    targetBSM.Add(WHBSM)
    othersignalBSM.Add(VBFBSM)
    
  VVHSM = VBFSM.Clone("VVHSM")
  VVHSM.Add(ZHSM)
  VVHSM.Add(WHSM)
  ffHSM = ggHSM.Clone("ffHSM")
  ffHSM.Add(ttHSM)
  ffHSM.Add(bbHSM)

  VVHBSM = VBFBSM.Clone("VVHBSM")
  VVHBSM.Add(ZHBSM)
  VVHBSM.Add(WHBSM)
  ffHBSM = ggHBSM.Clone("ffHBSM")
  ffHBSM.Add(ttHBSM)
  ffHBSM.Add(bbHBSM)

  ZZ = qqZZ.Clone("ZZ")
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

  targetSM.SetLineColor(4)
  targetSM.SetLineWidth(2)
  targetSM.SetFillStyle(0)
  othersignalSM.SetLineColor(ROOT.kOrange+10)
  othersignalSM.SetLineWidth(2)
  othersignalSM.SetFillStyle(0)

  targetBSM.SetLineColor(4)
  targetBSM.SetLineWidth(2)
  targetBSM.SetLineStyle(2)
  targetBSM.SetFillStyle(0)
  othersignalBSM.SetLineColor(ROOT.kOrange+10)
  othersignalBSM.SetLineWidth(2)
  othersignalBSM.SetLineStyle(2)
  othersignalBSM.SetFillStyle(0)

  histogramorderintegral = ffHSM, VVHSM, ffHBSM, VVHBSM, ZZ, ZX
  histogramorderstyle = othersignalSM, targetSM, othersignalBSM, targetBSM, ZZ, ZX
  normalization = {h: 0 for h in histogramorderintegral+(data,)}

  nfiles = 0
  g = glob.glob(os.path.join(config.plotsbasedir, "templateprojections/niceplots/fullrange", str(analysis), "*tagged", "D_bkg*.root"))
  for filename in g:
    nfiles += 1
    with TFile(filename) as f:
      for i, (hintegral, hstyle) in enumerate(itertools.izip(histogramorderintegral, histogramorderstyle)):
        proj = f.cprojections.GetListOfPrimitives()[1].GetHists()[i]
        try:
          assert proj.GetMarkerColor() == hstyle.GetMarkerColor(), (hstyle.GetMarkerColor(), proj.GetMarkerColor())
          assert proj.GetMarkerStyle() == hstyle.GetMarkerStyle(), (hstyle.GetMarkerStyle(), proj.GetMarkerStyle())
          assert proj.GetMarkerSize() == hstyle.GetMarkerSize(), (hstyle.GetMarkerSize(), proj.GetMarkerSize())
          assert proj.GetLineColor() == hstyle.GetLineColor(), (hstyle.GetLineColor(), proj.GetLineColor())
          assert proj.GetLineStyle() == hstyle.GetLineStyle(), (hstyle.GetLineStyle(), proj.GetLineStyle())
          assert proj.GetLineWidth() == hstyle.GetLineWidth(), (hstyle.GetLineWidth(), proj.GetLineWidth())
          assert proj.GetFillColor() == hstyle.GetFillColor(), (hstyle.GetFillColor(), proj.GetFillColor())
          assert proj.GetFillStyle() == hstyle.GetFillStyle(), (hstyle.GetFillStyle(), proj.GetFillStyle())
        except AssertionError:
          print i
          f.cprojections.GetListOfPrimitives()[1].ls()
          raise
        normalization[hintegral] += proj.Integral()

      datagraph = f.cprojections.GetListOfPrimitives()[6]
      normalization[data] += sum(y for _, y in itertools.izip(xrange(datagraph.GetN()), datagraph.GetY()))
  assert nfiles == 3, g

  normalization[ffHSM] -= normalization[VVHSM]
  normalization[ffHBSM] -= normalization[VVHBSM]
  normalization[VVHSM] -= normalization[ZZ]
  normalization[VVHBSM] -= normalization[ZZ]
  normalization[ZZ] -= normalization[ZX]

  for h, integral in normalization.iteritems():
    print h.h.GetName(), h.Integral(), integral
    if h == data:
      assert h.Integral() == integral, (h.Integral(), integral)
    elif h in (ffHSM, ffHBSM):
      assert abs(h.Integral() - integral) / (h.Integral() + integral) < 0.05, (h.Integral(), integral)
    elif h in (ZZ,):
      assert abs(h.Integral() - integral) / (h.Integral() + integral) < 1e-2, (h.Integral(), integral)
    else:
      assert abs(h.Integral() - integral) / (h.Integral() + integral) < 1e-4, (h.Integral(), integral)

  data = style.asymmerrorsfromhistogram(data.h, showemptyerrors=False)
  data.SetLineColor(1)
  data.SetMarkerColor(1)
  data.SetLineStyle(1)
  data.SetLineWidth(1)
  data.SetMarkerStyle(20)
  data.SetMarkerSize(1.2)

  ZZ.Add(ZX)
  targetSM.Add(ZZ)
  targetBSM.Add(ZZ)

  totalSM = othersignalSM.Clone("totalSM")
  totalSM.Add(targetSM)
  totalBSM = othersignalBSM.Clone("totalBSM")
  totalBSM.Add(targetBSM)

  assert abs(totalSM.Integral() / (ZZ.Integral() + VVHSM.Integral() + ffHSM.Integral()) - 1) < 1e-4, (totalSM.Integral(), ZZ.Integral() + VVHSM.Integral() + ffHSM.Integral())
  assert abs(totalSM.Integral() / (targetSM.Integral() + othersignalSM.Integral()) - 1) < 1e-4, (totalSM.Integral(), targetSM.Integral() + othersignalSM.Integral())
  assert abs(totalBSM.Integral() / (ZZ.Integral() + VVHBSM.Integral() + ffHBSM.Integral()) - 1) < 1e-4, (totalBSM.Integral(), ZZ.Integral() + VVHBSM.Integral() + ffHBSM.Integral())
  assert abs(totalBSM.Integral() / (targetBSM.Integral() + othersignalBSM.Integral()) - 1) < 1e-4, (totalBSM.Integral(), targetBSM.Integral() + othersignalBSM.Integral())

  hstack = ROOT.THStack()
  hstack.Add(ZZ.h, "hist")
  hstack.Add(ZX.h, "hist")
  hstack.Add(targetSM.h, "hist")
  hstack.Add(targetBSM.h, "hist")
  hstack.Add(totalSM.h, "hist")
  hstack.Add(totalBSM.h, "hist")

  ymax = style.ymax((hstack, "nostack"), (data, "P"))
  if VBForVH == "VBF": ymax = max(ymax, 12)
  hstack.SetMaximum(ymax)

  c = plotcopier.TCanvas("c", "", 8, 30, 800, 800)
  style.applycanvasstyle(c)
  c.SetBottomMargin(0.14)
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
  hstack.GetXaxis().SetTitleOffset(1.15)

  hstack.GetXaxis().CenterTitle()
  hstack.GetYaxis().CenterTitle()

  l = ROOT.TLegend(*(
    {
      "VBF": {
        Analysis("fa3"):   (.4, .55, .7, .9),
        Analysis("fa2"):   (.6, .55, .9, .9),
        Analysis("fL1"):   (.6, .55, .9, .9),
        Analysis("fL1Zg"): (.55, .55, .85, .9),
      }[analysis],
      "VH": (.6, .55, .9, .9),
    }[VBForVH]
  ))
  style.applylegendstyle(l)
  l.AddEntry(data, "Observed", "lp")
  l.AddEntry(totalSM.h, "Total SM", "l")
  l.AddEntry(targetSM.h, VBForVH+" SM", "f")
  l.AddEntry(totalBSM.h, "Total {} = 1".format(analysis.title()), "l")
  l.AddEntry(targetBSM.h, VBForVH+" {} = 1".format(analysis.title()), "f")
  l.AddEntry(ZZ.h, "ZZ/Z#gamma*", "f")
  l.AddEntry(ZX.h, "Z+X", "f")

  line = ROOT.TLine(0.5, 0, 0.5, ymax/2)
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
        for cut in "decay", "proddec", "bycategory", "m4l":
          if analysis != args.analysis is not None: continue
          if VBForVH != args.VBForVH is not None: continue
          if cut != args.cut is not None: continue
          publiccategorydiscriminants(analysis, VBForVH, pc, whichdbkg=cut)
