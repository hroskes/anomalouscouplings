from helperstuff import config
from helperstuff.enums import proddechypotheses, prodonlyhypotheses
from helperstuff.samples import ReweightingSample, Sample
import os
import ROOT
import helperstuff.style
ROOT.gStyle.SetCanvasDefW(678)
ROOT.gStyle.SetPadRightMargin(0.115)

fromsamples = [Sample("ZH", hypothesis, "160909") for hypothesis in prodonlyhypotheses]
tosamples = [ReweightingSample("ZH", hypothesis) for hypothesis in proddechypotheses]
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

for i, tosample in enumerate(tosamples):
    hstack = None
    del cache[:]
    hstack = ROOT.THStack(str(tosample).replace(" ", ""), "")
    legend = ROOT.TLegend(.6, .6, .9, .9)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    for fromsample in fromsamples:
        dirname = os.path.join(config.plotsbasedir, "weightplots/{}".format(str(fromsample).replace(" ", "")))
        try: os.makedirs(dirname)
        except: pass

        t = ROOT.TChain("candTree")
        t.Add(fromsample.withdiscriminantsfile())

        c1 = ROOT.TCanvas()
        c1.SetLogy()
        c1.SetLogx(False)
        hname = (str(fromsample)+str(tosample)).replace(" ", "")

        weightname = "reweightingweights[{}]".format(i)
        t.GetEntry(0)
        #weightname = tosample.weightname()
        #t.Draw("{}>>{}(200, 0, 45e-6)".format(weightname, hname))
        t.Draw("{}>>{}".format(weightname, hname))
        h = getattr(ROOT, hname)
        #if fromsample.hypothesis == "L1": continue
        hstack.Add(h)
        #h.SetLineColor(fromsample.color())
        h.SetDirectory(0)
        #cache.append(h)
        legend.AddEntry(h, str(fromsample.hypothesis), "l")
        for ext in "png", "eps", "root", "pdf":
            c1.SaveAs("{}/{}.{}".format(dirname, tosample.hypothesis, ext))
        c1.SetLogx(True)
        for ext in "png", "eps", "root", "pdf":
            c1.SaveAs("{}/{}_log.{}".format(dirname, tosample.hypothesis, ext))
    """
    print list(hstack.GetHists())
    print [h.Scale(1/h.Integral()) for h in hstack.GetHists()]
    hstack.Draw("hist nostack")
    legend.Draw()
    for ext in "png", "eps", "root", "pdf":
        c1.SaveAs("~/www/TEST/weightplots/{}.{}".format(tosample.hypothesis, ext))
    """
