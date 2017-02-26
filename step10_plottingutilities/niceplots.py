#!/usr/bin/env python
from helperstuff import config
from helperstuff import stylefunctions as style
from helperstuff.combinehelpers import Luminosity
from helperstuff.enums import Analysis, analyses, Channel, channels, EnumItem, MultiEnum, MyEnum, Production
from itertools import product
import os
from projections import EnrichStatus, enrichstatuses
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

def niceplots(*args, **kwargs):
    extradir = ""
    for kw, kwarg in kwargs.iteritems():
        if kw == "extradir":
            extradir = kwarg
        else:
            raise TypeError("Unknown kwarg {}={}".format(kw, kwarg))
    class NicePlots(MultiEnum):
        enums = (Analysis, Production, Category, EnrichStatus)
    info = NicePlots(*args)
    analysis, production, category, enrichstatus = info.analysis, info.production, info.category, info.enrichstatus
    dir = os.path.join(config.plotsbasedir, "templateprojections", "niceplots", extradir, enrichstatus.dirname(), str(analysis))

    previousplots = {(channel, production): os.path.join(config.plotsbasedir, "templateprojections", enrichstatus.dirname(), "rescalemixtures", "{}_{}".format(analysis, production), str(channel)) for channel, production in product(channels, productions)}

    for discriminant, title in ("Dbkg", "D_{bkg}"), (analysis.purediscriminant(), analysis.purediscriminant(True)), (analysis.mixdiscriminant(), analysis.mixdiscriminant(True)):

        ZX, ZZ, SM, BSM, mix, mixminus, data = HistWrapper(), HistWrapper(), HistWrapper(), HistWrapper(), HistWrapper(), HistWrapper(), HistWrapper()

        for channel, production in product(channels, productions):
            f = ROOT.TFile(os.path.join(previousplots[channel, production], discriminant+".root"))
            c = f.c1
            lst = c.GetListOfPrimitives()[1].GetHists()

            ZX += lst[6]
            ZZ += lst[4]
            ZZ += lst[5]
            SM += lst[0]
            BSM += lst[1]
            mix += lst[2]
            mixminus += lst[3]
            data += lst[7]

        ZX, ZZ, SM, BSM, mix, mixminus, data = ZX.h, ZZ.h, SM.h, BSM.h, mix.h, mixminus.h, data.h

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

        mixminus.SetLineColor(ROOT.kOrange+10)
        mixminus.SetLineStyle(2)
        mixminus.SetLineWidth(2)
        mixminus.Add(ZZ)  #which already has ZX

        data = style.asymmerrorsfromhistogram(data, showemptyerrors=False)
        data.SetLineColor(1)
        data.SetMarkerColor(1)
        data.SetLineStyle(1)
        data.SetLineWidth(1)
        data.SetMarkerStyle(20)
        data.SetMarkerSize(1.2)

        hstack = ROOT.THStack(discriminant, "")
        if discriminant == "Dbkg":
            l = ROOT.TLegend(0.23,0.57,0.61,0.90)
        else:
            l = ROOT.TLegend(0.20,0.57,0.58,0.90)
        #l.SetBorderSize(0)
        #l.SetFillStyle(0)
        style.applylegendstyle(l)

        hstack.Add(SM, "hist")
        if discriminant == analysis.mixdiscriminant() and analysis == "fa3" or discriminant in (analysis.mixdiscriminant(), analysis.purediscriminant()) and analysis == "fL1":
            hstack.Add(mix, "hist")
        elif discriminant in (analysis.mixdiscriminant(), analysis.purediscriminant()) and analysis == "fa2":
            hstack.Add(mixminus, "hist")
        else:
            hstack.Add(BSM, "hist")
        hstack.Add(ZZ, "hist")
        hstack.Add(ZX, "hist")

        l.AddEntry(data, "Observed", "ep")
        l.AddEntry(SM, "SM", "l")
        if discriminant == analysis.mixdiscriminant() and analysis == "fa3" or discriminant in (analysis.mixdiscriminant(), analysis.purediscriminant()) and analysis == "fL1":
            l.AddEntry(mix, analysis.title()+" = #plus 0.5", "l")
        elif discriminant in (analysis.mixdiscriminant(), analysis.purediscriminant()) and analysis == "fa2":
            l.AddEntry(mixminus, analysis.title()+" = #minus 0.5", "l")
        else:
            l.AddEntry(BSM, analysis.title()+" = 1", "l")
        l.AddEntry(ZZ, "ZZ/Z#gamma*", "f")
        l.AddEntry(ZX, "Z+X", "f")

        ymax = style.ymax((hstack, "nostack"), (data, "P"))
        hstack.SetMaximum(ymax*1.25)

        c = ROOT.TCanvas("c1", "", 8, 30, 800, 800)
        hstack.Draw("nostack")
        data.Draw("P")
        l.Draw()

        if discriminant == "D_CP_decay": hstack.GetXaxis().SetRangeUser(-.4, .4)

        style.applycanvasstyle(c)
        style.applyaxesstyle(hstack)
        style.cuttext(enrichstatus.cuttext())
        style.CMS("Preliminary", sum(float(Luminosity("fordata", production)) for production in productions))

        hstack.GetXaxis().SetTitle(title)
        hstack.GetYaxis().SetTitle(
                                   "Events / {:.2f}".format(
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
        for enrichstatus in enrichstatuses:
            p2015 = [p for p in config.productionsforcombine if p.year == 2015]
            p2016 = [p for p in config.productionsforcombine if p.year == 2016]
            assert len(p2015) == len(p2016) == 1
            niceplots(p2015, analysis, enrichstatus, extradir="2015")
            niceplots(p2016, analysis, enrichstatus, extradir="2016")
            niceplots(p2016, analysis, enrichstatus, extradir="")
            niceplots(config.productionsforcombine, analysis, enrichstatus, extradir="2015+2016")
