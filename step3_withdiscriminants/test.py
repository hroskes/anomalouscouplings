import ROOT
from helperstuff import style

f = ROOT.TFile("ggH0+160714.root")
t = f.candTree
c = ROOT.TCanvas()
t.Draw("p0plus_VAJHU*p0plus_m4l / (p0plus_VAJHU*p0plus_m4l + 1.34*bkg_VAMCFM*bkg_m4l)>>h", "(ZZMass > 105 && ZZMass < 140) * MC_weight_ggH_g1")
c.SaveAs("~/www/TEST/ggH_Dbkg.png")

h = ROOT.h
h.Scale(1/h.Integral())
print h.Integral(h.Fill(.5, 0), h.GetNbinsX())

f = ROOT.TFile("qqZZ160714.root")
t = f.candTree
c2 = ROOT.TCanvas()
t.Draw("p0plus_VAJHU*p0plus_m4l / (p0plus_VAJHU*p0plus_m4l + 1.34*bkg_VAMCFM*bkg_m4l)>>h2", "(ZZMass > 105 && ZZMass < 140) * MC_weight_qqZZ")
c2.SaveAs("~/www/TEST/qqZZ_Dbkg.png")

h2 = ROOT.h2
h2.Scale(1/h2.Integral())
print h2.Integral(h2.Fill(.5, 0), h2.GetNbinsX())




f = ROOT.TFile("ggH0+160714.root")
t = f.candTree
c = ROOT.TCanvas()
t.Draw("p0plus_VAJHU / (p0plus_VAJHU + 1.34*bkg_VAMCFM)>>h3", "(ZZMass > 105 && ZZMass < 140) * MC_weight_ggH_g1")
c.SaveAs("~/www/TEST/ggH_Dbkgkin.png")

h = ROOT.h3
h.Scale(1/h.Integral())
print h.Integral(h.Fill(.5, 0), h.GetNbinsX())

f = ROOT.TFile("qqZZ160714.root")
t = f.candTree
c2 = ROOT.TCanvas()
t.Draw("p0plus_VAJHU / (p0plus_VAJHU + 1.34*bkg_VAMCFM)>>h4", "(ZZMass > 105 && ZZMass < 140) * MC_weight_qqZZ")
c2.SaveAs("~/www/TEST/qqZZ_Dbkgkin.png")

h2 = ROOT.h4
h2.Scale(1/h2.Integral())
print h2.Integral(h2.Fill(.5, 0), h2.GetNbinsX())





f = ROOT.TFile("ggH0+160714.root")
t = f.candTree
c = ROOT.TCanvas()
t.Draw("p0plus_m4l / (p0plus_m4l + bkg_m4l)>>h5", "(ZZMass > 105 && ZZMass < 140) * MC_weight_ggH_g1")
c.SaveAs("~/www/TEST/ggH_Dbkgm4l.png")

h = ROOT.h5
h.Scale(1/h.Integral())
print h.Integral(h.Fill(.5, 0), h.GetNbinsX())

f = ROOT.TFile("qqZZ160714.root")
t = f.candTree
c2 = ROOT.TCanvas()
t.Draw("p0plus_m4l / (p0plus_m4l + bkg_m4l)>>h6", "(ZZMass > 105 && ZZMass < 140) * MC_weight_qqZZ")
c2.SaveAs("~/www/TEST/qqZZ_Dbkgm4l.png")

h2 = ROOT.h6
h2.Scale(1/h2.Integral())
print h2.Integral(h2.Fill(.5, 0), h2.GetNbinsX())





f = ROOT.TFile("ggH0+160714.root")
t = f.candTree
c = ROOT.TCanvas()
t.Draw("D_bkg_0plus>>h7", "(ZZMass > 105 && ZZMass < 140) * MC_weight_ggH_g1")
c.SaveAs("~/www/TEST/ggH_Dbkg_fixc.png")

h = ROOT.h7
h.Scale(1/h.Integral())
print h.Integral(h.Fill(.5, 0), h.GetNbinsX())

f = ROOT.TFile("qqZZ160714.root")
t = f.candTree
c2 = ROOT.TCanvas()
t.Draw("D_bkg_0plus>>h8", "(ZZMass > 105 && ZZMass < 140) * MC_weight_qqZZ")
c2.SaveAs("~/www/TEST/qqZZ_Dbkg_fixc.png")

h2 = ROOT.h8
h2.Scale(1/h2.Integral())
print h2.Integral(h2.Fill(.5, 0), h2.GetNbinsX())




