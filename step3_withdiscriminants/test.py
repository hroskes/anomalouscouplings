from helperstuff import config
from helperstuff import style
from helperstuff.enums import *
from helperstuff.samples import *
import ROOT

hstack = ROOT.THStack()
legend = ROOT.TLegend(.6, .5, .9, .9)
cache = []

ROOT.gErrorIgnoreLevel = ROOT.kError

c = ROOT.TCanvas()
for color, hypothesis in enumerate(prodonlyhypotheses, start=1):
    t = ROOT.TChain("candTree", "candTree")
    t.Add(Sample("VBF", hypothesis, "160729").withdiscriminantsfile())
    hname = "h{}".format(hypothesis)
    t.Draw("D_0minus_decay>>{}(100,-.5,.5)".format(hname), "MC_weight_VBF_g1", "hist")
    h = getattr(ROOT, hname)
    hstack.Add(h)
    cache.append(h)
    legend.AddEntry(h, str(hypothesis), "l")
    h.SetLineColor(color)
    print hypothesis, h.Integral()
    c.SaveAs("~/www/TEST/{}.png".format(hypothesis))

hstack.Draw("histnostack")
c.SaveAs("~/www/TEST/test.png")
