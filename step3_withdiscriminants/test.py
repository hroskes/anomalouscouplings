import ROOT
from helperstuff import style

f = ROOT.TFile("ggH0+160714.root")
t = f.candTree
c = ROOT.TCanvas()
t.Draw("D_bkg_0plus", "(ZZMass > 105 && ZZMass < 140) * MC_weight_ggH_g1")
c.SaveAs("~/www/TEST/test.png")

f = ROOT.TFile("ggZZ2e2mu160714.root")
t = f.candTree
c = ROOT.TCanvas()
t.Draw("D_bkg_0plus", "(ZZMass > 105 && ZZMass < 140) * MC_weight_ggZZ")
c.SaveAs("~/www/TEST/test2.png")
