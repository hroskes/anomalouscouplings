#!/usr/bin/env python

import ROOT
from helperstuff.CJLSTscripts import getDbkgkinConstant

f = ROOT.TFile("qqZZ161221.root")

t = f.candTree

h = ROOT.TH1F("h", "h", 50, 0, 1)

t.SetBranchStatus("*", 0)
t.SetBranchStatus("Z*Flav", 1)
t.SetBranchStatus("ZZMass", 1)
t.SetBranchStatus("p_GG_SIG_ghg2_1_ghz1_1_JHUGen", 1)
t.SetBranchStatus("p_QQB_BKG_MCFM", 1)
t.SetBranchStatus("MC_weight_qqZZ", 1)

for entry in t:
    ZZFlav = t.Z1Flav*t.Z2Flav
    ZZMass = t.ZZMass
    c = getDbkgkinConstant(ZZFlav, ZZMass)
    D = t.p_GG_SIG_ghg2_1_ghz1_1_JHUGen / (t.p_GG_SIG_ghg2_1_ghz1_1_JHUGen + t.p_QQB_BKG_MCFM*c)
    h.Fill(D, t.MC_weight_qqZZ)

c = ROOT.TCanvas()
h.Draw("hist")
c.SaveAs("~/www/TEST/test.png")
