from helperstuff.enums import hypotheses
import ROOT
from helperstuff.samples import Sample
import style
ROOT.gStyle.SetCanvasDefW(678)
ROOT.gStyle.SetPadRightMargin(0.115)

samples = [Sample("ggH", hypothesis) for hypothesis in hypotheses]
cache = []
"""
upperlimit = {
              "0+": 40,
              "0-": 40,
              "fa30.5": 40,
              "a2": 10,
              "fa20.5": 40,
              "L1": 40,
              "fL10.5": 25,
             }
"""

for i, tosample in enumerate(samples):
    hstack = None
    del cache[:]
    hstack = ROOT.THStack(str(tosample).replace(" ", ""), "")
    legend = ROOT.TLegend(.6, .6, .9, .9)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    for fromsample in samples:
        t = ROOT.TChain("candTree")
        t.Add(fromsample.withdiscriminantsfile())

        c1 = ROOT.TCanvas()
        c1.SetLogy()
        hname = (str(fromsample)+str(tosample)).replace(" ", "")

        weightname = "reweightingweights[{}]".format(i)
        t.GetEntry(0)
        #weightname = tosample.weightname()
        #t.Draw("{}>>{}(200, 0, 45e-6)".format(weightname, hname))
        t.Draw("{}>>{}".format(weightname, hname))
        h = getattr(ROOT, hname)
        #if fromsample.hypothesis == "L1": continue
        hstack.Add(h)
        h.SetLineColor(fromsample.color())
        h.SetDirectory(0)
        cache.append(h)
        legend.AddEntry(h, str(fromsample.hypothesis), "l")
    print list(hstack.GetHists())
    print [h.Scale(1/h.Integral()) for h in hstack.GetHists()]
    hstack.Draw("hist nostack")
    legend.Draw()
    for ext in "png", "eps", "root", "pdf":
        c1.SaveAs("~/www/TEST/weightplots/{}.{}".format(tosample.hypothesis, ext))
