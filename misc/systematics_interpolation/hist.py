#!/usr/bin/env python

import ROOT, rootoverloads

def g(h): return frozenset(x for x in xrange(1, h.GetNbinsX()+1) if h.GetBinContent(x))

f = ROOT.TFile("hzz4l_2e2muS_Untagged_2016.input.root")

hists = [
f.VH_g11g42g1prime21_negative,
f.VH_g11g42g1prime21_negative_CMS_scale_2e2muUp,
f.VH_g11g42g1prime21_negative_CMS_scale_2e2muDown,
f.VH_g11g42g1prime21_negative_CMS_res_2e2muUp,
f.VH_g11g42g1prime21_negative_CMS_res_2e2muDown,
]

sets = {g(_) for _ in hists}

#print frozenset.union(*sets) - frozenset.intersection(*sets)

#print [_.Integral() for _ in hists]

for h in f.VH_g11g42g1prime21_positive, f.VH_g11g42g1prime21_negative, f.VH_g11g42g1prime21_positive_CMS_scale_2e2muUp, f.VH_g11g42g1prime21_negative_CMS_scale_2e2muUp, f.VH_g11g42g1prime21_positive_CMS_scale_2e2muDown, f.VH_g11g42g1prime21_negative_CMS_scale_2e2muDown:
  print [h.GetBinContent(i) for i in xrange(1, h.GetNbinsX()+1)]
