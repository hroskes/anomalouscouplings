#!/usr/bin/env python
from helperstuff import config
from helperstuff import stylefunctions as style
from helperstuff.combinehelpers import Luminosity
from helperstuff.enums import Analysis, analyses, categories, Category, Channel, channels, EnumItem, MultiEnum, MyEnum, Production
from helperstuff.templates import TemplatesFile
from itertools import product
import os
from projections import EnrichStatus, enrichstatuses
import ROOT

assert len(config.productionsforcombine) == 1
production = config.productionsforcombine[0]

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
    dir = os.path.join(config.plotsbasedir, "templateprojections", "niceplots", extradir, enrichstatus.dirname(), str(analysis), str(category))

    previousplots = {channel: os.path.join(config.plotsbasedir, "templateprojections", enrichstatus.dirname(), "{}_{}".format(analysis, production), str(category), str(channel)) for channel in channels}

    tf = TemplatesFile("2e2mu", "ggh", category, production, analysis)
    for discriminant, title, _, _, _ in tf.discriminants:

        ZX, ZZ, SM, BSM, mix, mixminus, data = HistWrapper(), HistWrapper(), HistWrapper(), HistWrapper(), HistWrapper(), HistWrapper(), HistWrapper()

        for channel in channels:
            f = ROOT.TFile(os.path.join(previousplots[channel], discriminant+".root"))
            c = f.c1
            lst = c.GetListOfPrimitives()[1].GetHists()

            if category == "VBFtagged" or category == "VHHadrtagged":
                for h in lst:
                    h.Rebin(5)

            ZX += lst[7]
            ZZ += lst[4]
            ZZ += lst[5]
            ZZ += lst[6]
            SM += lst[0]
            BSM += lst[1]
            mix += lst[2]
            mixminus += lst[3]
            try:
                data += lst[8]
            except IndexError:
                return

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

        l.AddEntry(data, "Observed", "ep")
        hstack.Add(SM, "hist")
        l.AddEntry(SM, "SM", "l")
        if discriminant == tf.mixdiscriminant and analysis == "fa3" or discriminant in (tf.mixdiscriminant, tf.purediscriminant) and analysis == "fL1":
            hstack.Add(mix, "hist")
            l.AddEntry(mix, analysis.title()+" = #plus 0.5", "l")
        elif discriminant in (tf.mixdiscriminant, tf.purediscriminant) and analysis == "fa2":
            hstack.Add(mixminus, "hist")
            l.AddEntry(mixminus, analysis.title()+" = #minus 0.5", "l")
        else:
            hstack.Add(BSM, "hist")
            l.AddEntry(BSM, analysis.title()+" = 1", "l")
        hstack.Add(ZZ, "hist")
        l.AddEntry(ZZ, "ZZ/Z#gamma*", "f")
        hstack.Add(ZX, "hist")
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
        style.CMS("Preliminary", float(Luminosity("fordata", production)))

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
            for category in categories:
                niceplots(production, category, analysis, enrichstatus)
