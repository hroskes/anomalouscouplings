from helperstuff import config
from helperstuff import stylefunctions as style
from helperstuff.enums import Analysis, analyses, Channel, channels, MultiEnum, Production
import os
import ROOT

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

    if config.unblinddata:
        blinddir = "fullrange"
    else:
        blinddir = "blind"
    previousplots = {channel: os.path.join(config.plotsbasedir, "templateprojections", blinddir, "rescalemixtures", "{}_{}".format(analysis, production), str(channel)) for channel in channels}

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
        ZX.SetLineWidth(2)
        ZX.SetFillStyle(1001)

        ZZ.SetLineColor(1)
        ZZ.SetFillColor(ROOT.kAzure-9)
        ZZ.SetLineWidth(2)
        ZX.SetFillStyle(1001)
        ZZ.Add(ZX)

        SM.SetLineColor(ROOT.kOrange+10)
        SM.SetLineWidth(2)
        SM.Add(ZZ)  #which already has ZX

        BSM.SetLineColor(ROOT.kOrange+10)
        BSM.SetLineStyle(2)
        BSM.SetLineWidth(2)
        BSM.Add(ZZ)  #which already has ZX

        mix.SetLineColor(ROOT.kOrange+10)
        mix.SetLineStyle(2)
        mix.SetLineWidth(2)
        mix.Add(ZZ)  #which already has ZX

        data = style.asymmerrorsfromhistogram(data, showemptyerrors=True)
        data.SetLineColor(1)
        data.SetMarkerColor(1)
        data.SetLineStyle(1)
        data.SetLineWidth(1)
        data.SetMarkerStyle(20)
        data.SetMarkerSize(1.2)

        hstack = ROOT.THStack(discriminant, "")
        l = ROOT.TLegend(.2, .5, .4, .95)
        l.SetBorderSize(0)
        l.SetFillStyle(0)

        hstack.Add(SM, "hist")
        if discriminant == analysis.mixdiscriminant():
            hstack.Add(mix, "hist")
        else:
            hstack.Add(BSM, "hist")
        hstack.Add(ZZ, "hist")
        hstack.Add(ZX, "hist")

        l.AddEntry(data, "data", "ep")
        l.AddEntry(SM, "SM", "l")
        if discriminant == analysis.mixdiscriminant():
            l.AddEntry(mix, analysis.title()+"=0.5", "l")
        else:
            l.AddEntry(BSM, analysis.title()+"=1", "l")
        l.AddEntry(ZZ, "q#bar{q}/gg#rightarrowZZ", "f")
        l.AddEntry(ZX, "Z+X", "f")

        ymax = style.ymax((hstack, "nostack"), (data, "P"))
        hstack.SetMaximum(ymax*1.25)

        c = ROOT.TCanvas("c1", "", 8, 30, 800, 800)
        hstack.Draw("nostack")
        data.Draw("P")
        l.Draw()

        style.applycanvasstyle(c)
        style.applyaxesstyle(hstack)

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
