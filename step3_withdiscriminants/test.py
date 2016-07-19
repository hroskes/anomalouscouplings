import ROOT
import style

ROOT.gROOT.LoadMacro("test.C")
from math import *
def getDbkgkinConstant(ZZflav, ZZMass):
  par=[
    0.775,
    -0.565,
    70.,
    5.90,
    -0.235,
    130.1,
    13.25,
    -0.33,
    191.04,
    16.05,
    187.47,
    -0.58,
    1700.,
    400.,
  ]
  if abs(ZZflav)==121*121: par[11]=-0.68
  kappa = (
    par[0]
    +par[1]*exp(-pow(((ZZMass-par[2])/par[3]), 2))
    +par[4]*exp(-pow(((ZZMass-par[5])/par[6]), 2))
    +par[7]*(
    exp(-pow(((ZZMass-par[8])/par[9]), 2))*(ZZMass<par[8])
    + exp(-pow(((ZZMass-par[8])/par[10]), 2))*(ZZMass>=par[8])
    +par[11]*exp(-pow(((ZZMass-par[12])/par[13]), 2))
    )
  )
  constant = kappa/(1.-kappa)
  return constant


t = ROOT.TChain("ZZTree/candTree")
print t
t.Add("root://lxcms03//data3/Higgs/160714/ggH1500/ZZ4lAnalysis.root")
print t
h = ROOT.TH1F("ggH", "ggH", 100, 0, 1)
lenght = t.GetEntries()
for i, entry in enumerate(t):
    if 1200 < t.ZZMass < 2000 and t.Z1Flav*t.Z2Flav == 11**4:
        D_bkg_kin = t.p0plus_VAJHU / (t.p0plus_VAJHU + getDbkgkinConstant(11**4, t.ZZMass))
        h.Fill(D_bkg_kin)
    if i % 1000 == 0: print i, "/", lenght

t = ROOT.TChain("ZZTree/candTree")
print t
t.Add("root://lxcms03//data3/Higgs/160714/ZZTo4l/ZZ4lAnalysis.root")
print t
h2 = ROOT.TH1F("qqZZ", "qqZZ", 100, 0, 1)
h2.SetLineColor(2)
lenght = t.GetEntries()
n=0
for i, entry in enumerate(t):
    if 1200 < t.ZZMass < 2000 and t.Z1Flav*t.Z2Flav == 11**4:
        D_bkg_kin = t.p0plus_VAJHU / (t.p0plus_VAJHU + getDbkgkinConstant(11**4, t.ZZMass))
        h2.Fill(D_bkg_kin)
        n+=1
    if i % 1000 == 0: print i, "/", lenght

h.Scale(1/h.Integral())
h2.Scale(1/h2.Integral())

hstack = ROOT.THStack("hs", "hs")
hstack.Add(h)
hstack.Add(h2)
c1 = ROOT.TCanvas()
hstack.Draw("nostack")
c1.SaveAs("~/www/TEST/test.png")
print n
