from helperstuff import config
from helperstuff import style
from helperstuff.enums import *
from helperstuff.samples import *
import ROOT
import os

hstack = ROOT.THStack()
legend = ROOT.TLegend(.6, .5, .9, .9)
cache = []

ROOT.gErrorIgnoreLevel = ROOT.kError

c = ROOT.TCanvas()
for color, hypothesis in enumerate(["L1", "fL1prod0.5", "0-", "fa3prod0.5", "a2", "fa2prod0.5", "0+"], start=1):
    t = ROOT.TChain("candTree", "candTree")
    t.Add(Sample("ZH", hypothesis, "160909").withdiscriminantsfile())
    hname = "h{}".format(hypothesis)
    #weight = "MC_weight_ZH_g1g4_prod * (D_CP_ZH>-998)"
    weight = "(D_CP_ZH_hadronic>-998)"
    t.Draw("D_g13_g21_ZHdecay_hadronic_prime>>{}(50,-1,1)".format(hname), weight, "hist")
#    t.Draw("D_g12_g1prime22_ZHdecay:D_g1g1prime2_ZHdecay>>{}".format(hname), weight, "SCAT")
#    t.Draw("D_CP_decay>>{}(50,-.5,.5)".format(hname), weight, "hist")
#    t.Draw("D_CP_ZH>>{}(50,-1,1)".format(hname), weight+"*(D_2jet_0plus>-1)", "hist")
    h = getattr(ROOT, hname)
    if isinstance(h, ROOT.TH1) and not isinstance(h, ROOT.TH2):
      h.SetBinContent(h.GetNbinsX(), h.GetBinContent(h.GetNbinsX()+1) + h.GetBinContent(h.GetNbinsX()))
      h.SetBinContent(1, h.GetBinContent(0) + h.GetBinContent(1))
    hstack.Add(h)
    cache.append(h)
    legend.AddEntry(h, str(hypothesis), "l")
    h.SetLineColor(color)
    h.SetMarkerStyle(1)
    print hypothesis, h.Integral()
    try:
        os.makedirs(os.path.join(config.plotsbasedir, "TEST"))
    except OSError:
        pass
    c.SaveAs(os.path.join(config.plotsbasedir, "TEST", "{}.png".format(hypothesis)))

hstack.Draw("histnostack")
c.SaveAs(os.path.join(config.plotsbasedir, "TEST", "test.png"))
