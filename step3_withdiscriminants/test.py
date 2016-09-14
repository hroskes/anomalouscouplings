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
    weight = "MC_weight_ZH_g1g4_prod"
    t.Draw("D_CP_ZH_hadronic>>{}(50,-.1,.1)".format(hname), weight, "hist")
#    t.Draw("D_g1prime2_VBF:D_g1g1prime2_VBF>>{}".format(hname), weight, "SCAT")
#    t.Draw("D_CP_decay>>{}(50,-.5,.5)".format(hname), weight, "hist")
#    t.Draw("D_CP_VBF>>{}(50,-1,1)".format(hname), weight+"*(D_2jet_0plus>-1)", "hist")
    h = getattr(ROOT, hname)
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
