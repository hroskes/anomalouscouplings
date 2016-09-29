from helperstuff import config
from helperstuff import style
from helperstuff.discriminants import discriminant
from helperstuff.enums import *
from helperstuff.samples import *
import ROOT
import os

#========================
#inputs
#weight, bins, min, max can be None
productionmode = "WH"
disc           = "D_0minus_WH_hadronic"
weight         = "MC_weight_{}_g1".format(productionmode)
bins           = None
min            = None
max            = None
#========================


hstack = ROOT.THStack()
legend = ROOT.TLegend(.6, .5, .9, .9)
cache = []

ROOT.gErrorIgnoreLevel = ROOT.kError

discname, title, discbins, discmin, discmax = discriminant(disc)
if bins is None:
    bins = discbins
if min is None:
    min = discmin
if max is None:
    max = discmax

c = ROOT.TCanvas()
hs = {}
for color, hypothesis in enumerate(["L1", "fL1prod0.5", "0-", "fa3prod0.5", "a2", "fa2prod0.5", "0+"], start=1):
#    if hypothesis in ["0-", "fa3prod0.5"]: continue
    t = ROOT.TChain("candTree", "candTree")
    sample = Sample(productionmode, hypothesis, "160928")
    t.Add(sample.withdiscriminantsfile())
    hname = "h{}".format(hypothesis)

    weightname = weight if weight is not None else sample.weightname()

    wt = "({}>-998)*{}".format(discname, weightname)
    t.Draw("{}>>{}({},{},{})".format(discname, hname, bins, min, max), wt, "hist")
    h = hs[hypothesis] = getattr(ROOT, hname)
    if isinstance(h, ROOT.TH1) and not isinstance(h, ROOT.TH2):
      h.SetBinContent(h.GetNbinsX(), h.GetBinContent(h.GetNbinsX()+1) + h.GetBinContent(h.GetNbinsX()))
      h.SetBinContent(1, h.GetBinContent(0) + h.GetBinContent(1))
    hstack.Add(h)
    cache.append(h)
    legend.AddEntry(h, str(hypothesis), "l")
    h.SetLineColor(color)
    h.SetMarkerStyle(1)
    print "{:10} {:.3g}".format(hypothesis, h.Integral())
    try:
        os.makedirs(os.path.join(config.plotsbasedir, "TEST"))
    except OSError:
        pass
    c.SaveAs(os.path.join(config.plotsbasedir, "TEST", "{}.png".format(hypothesis)))

hstack.Draw("histnostack")
c.SaveAs(os.path.join(config.plotsbasedir, "TEST", "test.png"))

"""
hint = hs["0+"].Clone("hint")
hint.Add(hs["0-"])
hint.Scale(-.5)
hint.Add(hs["fa3prod0.5"])

hint.Draw("hist")
#c.SaveAs(os.path.join(config.plotsbasedir, "TEST", "a1a3int.png"))
"""
