from helperstuff import config
from helperstuff.enums import Analysis, analyses, Channel, channels, MultiEnum, Production
import os
import ROOT
import style

class HistWrapper(object):
    def __init__(self):
        self.h = 0
    def __add__(self, other):
        if self.h == 0:
            self.h = other.Clone()
            self.h.SetDirectory(0)
        else:
            self.h.Add(other)
        return self

def niceplots(*args):
    class NicePlots(MultiEnum):
        enums = (Analysis, Production)
    info = NicePlots(*args)
    analysis, production = info.analysis, info.production
    dir = os.path.join(config.plotsbasedir, "templateprojections", "niceplots", str(analysis))
    previousplots = {channel: os.path.join(config.plotsbasedir, "templateprojections", "blind", "rescalemixtures", "{}_{}".format(analysis, production), str(channel)) for channel in channels}

    for discriminant, title in ("Dbkg", "D_{bkg}"), (analysis.purediscriminant(), analysis.purediscriminant(True)), (analysis.mixdiscriminant(), analysis.mixdiscriminant(True)):

        ZX, ZZ, SM, BSM, mix, data = HistWrapper(), HistWrapper(), HistWrapper(), HistWrapper(), HistWrapper(), HistWrapper()

        for channel in channels:
            f = ROOT.TFile(os.path.join(previousplots[channel], discriminant+".root"))
            c = f.c1
            lst = c.GetListOfPrimitives()[1].GetHists()

            ZX += lst[6]
            ZZ += lst[4]
            ZZ += lst[5]
            SM += lst[0]
            BSM += lst[1]
            mix += lst[2]
            data += lst[7]

        ZX, ZZ, SM, BSM, mix, data = ZX.h, ZZ.h, SM.h, BSM.h, mix.h, data.h

        ZX.SetLineColor(1)
        ZX.SetFillColor(ROOT.TColor.GetColor("#669966"))
        ZX.SetLineWidth(3)

        ZZ.SetLineColor(1)
        ZZ.SetFillColor(ROOT.kAzure-9)
        ZZ.SetLineWidth(3)
        ZZ.Add(ZX)

        SM.SetLineColor(2)
        SM.SetLineWidth(3)
        SM.Add(ZZ)  #which already has ZX

        BSM.SetLineColor(2)
        BSM.SetLineStyle(2)
        BSM.SetLineWidth(3)
        BSM.Add(ZZ)  #which already has ZX

        mix.SetLineColor(ROOT.kViolet)
        mix.SetLineStyle(3)
        mix.SetLineWidth(3)
        mix.Add(ZZ)  #which already has ZX

        data.SetLineColor(1)
        data.SetMarkerColor(1)
        data.SetLineStyle(1)
        data.SetLineWidth(3)
        data.SetMarkerStyle(20)

        hstack = ROOT.THStack(discriminant, "")
        l = ROOT.TLegend(.2, .5, .4, .95)
        l.SetBorderSize(0)
        l.SetFillStyle(0)

        for h in ZZ, ZX, SM, BSM, mix:
            hstack.Add(h, "hist")
        hstack.Add(data, "P")
        l.AddEntry(data, "data", "lp")
        l.AddEntry(SM, "SM", "l"),
        l.AddEntry(BSM, analysis.title()+"=1", "l"),
        l.AddEntry(mix, analysis.title()+"=0.5", "l"),
        l.AddEntry(ZZ, "q#bar{q}/gg#rightarrowZZ", "f"),
        l.AddEntry(ZX, "Z+X", "f"),

        c = ROOT.TCanvas()
        hstack.Draw("nostack")
        l.Draw()

        hstack.GetXaxis().SetTitle(title)
        hstack.GetYaxis().SetTitle(
                                   "number of events / {}".format(
                                                                  (hstack.GetXaxis().GetXmax() - hstack.GetXaxis().GetXmin()) / hstack.GetXaxis().GetNbins()
                                                                 )
                                  )

        try:
            os.makedirs(dir)
        except OSError:
            pass
        for ext in "png eps root pdf".split():
            c.SaveAs(os.path.join(dir, discriminant+"."+ext))

if __name__ == "__main__":
    for analysis in analyses:
        niceplots(analysis, "160714")
