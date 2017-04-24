#!/usr/bin/env python
import ROOT, helperstuff.style

from helperstuff.combinehelpers import gettemplate

t = [
     gettemplate("ggH", "SM", "170222", "fa3", "VBFtagged", "2e2mu"),
     gettemplate("ggH", "0-", "170222", "fa3", "VBFtagged", "2e2mu"),
     gettemplate("ggH", "g11gi1", "170222", "fa3", "VBFtagged", "2e2mu").Clone("plus"),
     gettemplate("ggH", "g11gi1", "170222", "fa3", "VBFtagged", "2e2mu").Clone("minus"),
    ]

t[2].Add(t[0])
t[2].Add(t[1])

t[3].Scale(-1)
t[3].Add(t[0])
t[3].Add(t[1])

hstack = ROOT.THStack("hs", "hs")
cache = []
c = ROOT.TCanvas()
for color, _ in enumerate(t, start=1):
    _.Scale(1/_.Integral())
    proj = _.ProjectionY()
    proj.SetLineColor(color)
    hstack.Add(proj)
    cache.append(proj)
hstack.Draw("hist nostack")
c.SaveAs("~/www/TEST/test.png")

for i in range(1, t[0].GetNbinsY()+1):
    print i, cache[3].GetBinContent(i) - cache[3].GetBinContent(t[0].GetNbinsY()+1-i)
