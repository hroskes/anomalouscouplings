#!/usr/bin/env python

import os

import ROOT

from helperstuff import config, style
from helperstuff.discriminants import discriminant
from helperstuff.samples import Sample
from helperstuff.utilities import mkdir_p, PlotCopier, TFile

def getsample(productionmode, hypothesis):
  args = "190821_2017",
  if productionmode == "VBF" and hypothesis == "0-": args = "190821_2017", "ext1"
  if productionmode == "VBF" and hypothesis == "a2": args = "190821_2017", "ext1"
  if productionmode == "VBF" and hypothesis == "L1": args = "190821_2018",
  if productionmode == "VBF" and hypothesis == "L1Zg": args = "190821_2018",
  if productionmode == "ZH" and hypothesis == "L1Zg": args = "190821_2018",
  return Sample(productionmode, hypothesis, *args)

def plot(discnominal, discsystematic, productionmode, cutnominal, cutsystematic, plotname, plotcopier=ROOT):
  c = plotcopier.TCanvas()

  name, _, bins, min, max, _, formula = discriminant(discnominal)
  systname, _, systbins, systmin, systmax, _, systformula = discriminant(discsystematic)
  assert bins == systbins
  assert min == systmin
  assert max == systmax

  hists = []

  for color, hypothesis in enumerate(("0+", "0-", "a2", "L1", "L1Zg"), start=1):
    if color == 5: color = 6
    filename = getsample(productionmode, hypothesis).withdiscriminantsfile()
    with TFile(filename) as f:
      t = f.candTree
      if not t.GetEntries(): raise ValueError("no entries " + filename)
      t.Draw("{}>>h{}({},{},{})".format(formula, color, bins, min, max), cutnominal)
      t.Draw("{}>>h{}syst({},{},{})".format(systformula, color, bins, min, max), cutsystematic)

      hsyst = ROOT.gDirectory.Get("h{}syst".format(color))
      hsyst.Sumw2()
      h = ROOT.gDirectory.Get("h{}".format(color))
      h.Sumw2()

      ratio = hsyst
      ratio.Divide(h)

      ratio.SetLineColor(color)
      ratio.SetMarkerColor(color)
      ratio.SetMarkerStyle(21)
      ratio.SetMarkerSize(.7)
      ratio.SetDirectory(0)

      hists.append(ratio)

  hstack = ROOT.THStack("hstack", "hstack")

  for h in hists: hstack.Add(h)
  hstack.Draw("nostack PE")
  hstack.SetMinimum(hstack.GetMinimum("nostack PE"))

  saveasdir = os.path.join(config.plotsbasedir, "xchecks", "ACsystematics")
  mkdir_p(saveasdir)
  for ext in "png", "pdf", "root", "C":
    c.SaveAs(os.path.join(saveasdir, plotname+"."+ext))

if __name__ == "__main__":
  with PlotCopier() as pc:
    plot("D_0minus_VBFdecay_3bins", "D_0minus_VBFdecay_3bins_JECUp", "VBF", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==2", "category_0P_or_0M_or_a2_or_L1_or_L1Zg_JECUp==2", "D_0minus_VBFdecay_JECUp", pc)
    plot("D_0minus_VBFdecay_3bins", "D_0minus_VBFdecay_3bins_JECDn", "VBF", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==2", "category_0P_or_0M_or_a2_or_L1_or_L1Zg_JECDn==2", "D_0minus_VBFdecay_JECDn", pc)
    plot("D_0hplus_VBFdecay_3bins", "D_0hplus_VBFdecay_3bins_JECUp", "VBF", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==2", "category_0P_or_0M_or_a2_or_L1_or_L1Zg_JECUp==2", "D_0hplus_VBFdecay_JECUp", pc)
    plot("D_0hplus_VBFdecay_3bins", "D_0hplus_VBFdecay_3bins_JECDn", "VBF", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==2", "category_0P_or_0M_or_a2_or_L1_or_L1Zg_JECDn==2", "D_0hplus_VBFdecay_JECDn", pc)
    plot("D_L1_VBFdecay_3bins", "D_L1_VBFdecay_3bins_JECUp", "VBF", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==2", "category_0P_or_0M_or_a2_or_L1_or_L1Zg_JECUp==2", "D_L1_VBFdecay_JECUp", pc)
    plot("D_L1_VBFdecay_3bins", "D_L1_VBFdecay_3bins_JECDn", "VBF", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==2", "category_0P_or_0M_or_a2_or_L1_or_L1Zg_JECDn==2", "D_L1_VBFdecay_JECDn", pc)
    plot("D_L1Zg_VBFdecay_3bins", "D_L1Zg_VBFdecay_3bins_JECUp", "VBF", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==2", "category_0P_or_0M_or_a2_or_L1_or_L1Zg_JECUp==2", "D_L1Zg_VBFdecay_JECUp", pc)
    plot("D_L1Zg_VBFdecay_3bins", "D_L1Zg_VBFdecay_3bins_JECDn", "VBF", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==2", "category_0P_or_0M_or_a2_or_L1_or_L1Zg_JECDn==2", "D_L1Zg_VBFdecay_JECDn", pc)
    plot("D_int_VBF_new_2bins", "D_int_VBF_new_2bins_JECUp", "VBF", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==2", "category_0P_or_0M_or_a2_or_L1_or_L1Zg_JECUp==2", "D_int_VBF_JECUp", pc)
    plot("D_int_VBF_new_2bins", "D_int_VBF_new_2bins_JECDn", "VBF", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==2", "category_0P_or_0M_or_a2_or_L1_or_L1Zg_JECDn==2", "D_int_VBF_JECDn", pc)
    plot("D_bkg_VBFdecay_3bins", "D_bkg_VBFdecay_3bins_JECUp", "VBF", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==2", "category_0P_or_0M_or_a2_or_L1_or_L1Zg_JECUp==2", "D_bkg_VBFdecay_JECUp", pc)
    plot("D_bkg_VBFdecay_3bins", "D_bkg_VBFdecay_3bins_JECDn", "VBF", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==2", "category_0P_or_0M_or_a2_or_L1_or_L1Zg_JECDn==2", "D_bkg_VBFdecay_JECDn", pc)

    plot("D_0minus_HadVHdecay_3bins", "D_0minus_HadVHdecay_3bins_JECUp", "ZH", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==4", "category_0P_or_0M_or_a2_or_L1_or_L1Zg_JECUp==4", "D_0minus_HadVHdecay_JECUp", pc)
    plot("D_0minus_HadVHdecay_3bins", "D_0minus_HadVHdecay_3bins_JECDn", "ZH", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==4", "category_0P_or_0M_or_a2_or_L1_or_L1Zg_JECDn==4", "D_0minus_HadVHdecay_JECDn", pc)
    plot("D_0hplus_HadVHdecay_3bins", "D_0hplus_HadVHdecay_3bins_JECUp", "ZH", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==4", "category_0P_or_0M_or_a2_or_L1_or_L1Zg_JECUp==4", "D_0hplus_HadVHdecay_JECUp", pc)
    plot("D_0hplus_HadVHdecay_3bins", "D_0hplus_HadVHdecay_3bins_JECDn", "ZH", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==4", "category_0P_or_0M_or_a2_or_L1_or_L1Zg_JECDn==4", "D_0hplus_HadVHdecay_JECDn", pc)
    plot("D_L1_HadVHdecay_3bins", "D_L1_HadVHdecay_3bins_JECUp", "ZH", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==4", "category_0P_or_0M_or_a2_or_L1_or_L1Zg_JECUp==4", "D_L1_HadVHdecay_JECUp", pc)
    plot("D_L1_HadVHdecay_3bins", "D_L1_HadVHdecay_3bins_JECDn", "ZH", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==4", "category_0P_or_0M_or_a2_or_L1_or_L1Zg_JECDn==4", "D_L1_HadVHdecay_JECDn", pc)
    plot("D_L1Zg_HadVHdecay_3bins", "D_L1Zg_HadVHdecay_3bins_JECUp", "ZH", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==4", "category_0P_or_0M_or_a2_or_L1_or_L1Zg_JECUp==4", "D_L1Zg_HadVHdecay_JECUp", pc)
    plot("D_L1Zg_HadVHdecay_3bins", "D_L1Zg_HadVHdecay_3bins_JECDn", "ZH", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==4", "category_0P_or_0M_or_a2_or_L1_or_L1Zg_JECDn==4", "D_L1Zg_HadVHdecay_JECDn", pc)
    plot("D_int_HadVH_new_2bins", "D_int_HadVH_new_2bins_JECUp", "ZH", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==4", "category_0P_or_0M_or_a2_or_L1_or_L1Zg_JECUp==4", "D_int_HadVH_JECUp", pc)
    plot("D_int_HadVH_new_2bins", "D_int_HadVH_new_2bins_JECDn", "ZH", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==4", "category_0P_or_0M_or_a2_or_L1_or_L1Zg_JECDn==4", "D_int_HadVH_JECDn", pc)
    plot("D_bkg_HadVHdecay_3bins", "D_bkg_HadVHdecay_3bins_JECUp", "ZH", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==4", "category_0P_or_0M_or_a2_or_L1_or_L1Zg_JECUp==4", "D_bkg_HadVHdecay_JECUp", pc)
    plot("D_bkg_HadVHdecay_3bins", "D_bkg_HadVHdecay_3bins_JECDn", "ZH", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==4", "category_0P_or_0M_or_a2_or_L1_or_L1Zg_JECDn==4", "D_bkg_HadVHdecay_JECDn", pc)

    plot("D_bkg_3bins", "D_bkg_ScaleUp_3bins", "VBF", "1", "1", "D_bkg_ScaleUp_3bins", pc)
    plot("D_bkg_3bins", "D_bkg_ScaleDown_3bins", "VBF", "1", "1", "D_bkg_ScaleDown_3bins", pc)
    plot("D_bkg_3bins", "D_bkg_ResUp_3bins", "VBF", "1", "1", "D_bkg_ResUp_3bins", pc)
    plot("D_bkg_3bins", "D_bkg_ResDown_3bins", "VBF", "1", "1", "D_bkg_ResDown_3bins", pc)

    plot("D_bkg_VBFdecay_3bins", "D_bkg_VBFdecay_ScaleUp_3bins", "VBF", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==2", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==2", "D_bkg_VBFdecay_ScaleUp_3bins", pc)
    plot("D_bkg_VBFdecay_3bins", "D_bkg_VBFdecay_ScaleDown_3bins", "VBF", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==2", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==2", "D_bkg_VBFdecay_ScaleDown_3bins", pc)
    plot("D_bkg_VBFdecay_3bins", "D_bkg_VBFdecay_ResUp_3bins", "VBF", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==2", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==2", "D_bkg_VBFdecay_ResUp_3bins", pc)
    plot("D_bkg_VBFdecay_3bins", "D_bkg_VBFdecay_ResDown_3bins", "VBF", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==2", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==2", "D_bkg_VBFdecay_ResDown_3bins", pc)

    plot("D_bkg_HadVHdecay_3bins", "D_bkg_HadVHdecay_ScaleUp_3bins", "ZH", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==4", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==4", "D_bkg_HadVHdecay_ScaleUp_3bins", pc)
    plot("D_bkg_HadVHdecay_3bins", "D_bkg_HadVHdecay_ScaleDown_3bins", "ZH", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==4", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==4", "D_bkg_HadVHdecay_ScaleDown_3bins", pc)
    plot("D_bkg_HadVHdecay_3bins", "D_bkg_HadVHdecay_ResUp_3bins", "ZH", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==4", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==4", "D_bkg_HadVHdecay_ResUp_3bins", pc)
    plot("D_bkg_HadVHdecay_3bins", "D_bkg_HadVHdecay_ResDown_3bins", "ZH", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==4", "category_0P_or_0M_or_a2_or_L1_or_L1Zg==4", "D_bkg_HadVHdecay_ResDown_3bins", pc)
