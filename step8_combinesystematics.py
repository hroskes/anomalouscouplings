#!/usr/bin/env python
"""
combine ScaleUp, ResUp --> ScaleResUp
ZX systematics templates from qqZZ
"""
from itertools import izip

import ROOT

import helperstuff.rootoverloads.histogramfloor

from helperstuff import config

from helperstuff.combinehelpers import gettemplate
from helperstuff.enums import analyses, categories, channels, productions
from helperstuff.templates import IntTemplate, Template, TemplatesFile
from helperstuff.utilities import tfiles

def combinesystematics(channel, analysis, production, category):
    thetemplatesfiles = []

    if config.applym4lshapesystematics:
        ScaleAndRes = {}
        for _ in "ggh", "vbf", "zh", "wh", "tth":
            for syst in "ScaleUp", "ScaleDown", "ResUp", "ResDown":
                if _ == "ggh" and category == "Untagged": continue
                ScaleAndRes[_,syst] = TemplatesFile(channel, _, syst, analysis, production, category)
        thetemplatesfiles += ScaleAndRes.values()
        if config.combinem4lshapesystematics:
            ScaleResUp, ScaleResDown = {}, {}
            for _ in "ggh", "vbf", "zh", "wh", "tth":
                ScaleResUp[_] = TemplatesFile(channel, _, "ScaleResUp", analysis, production, category)
                ScaleResDown[_] = TemplatesFile(channel, _, "ScaleResDown", analysis, production, category)
            thetemplatesfiles += ScaleResUp.values() + ScaleResDown.values()

    if config.applyZXshapesystematics:
        ZXUp = TemplatesFile(channel, "bkg", "ZXUp", analysis, production, category)
        ZXDown = TemplatesFile(channel, "bkg", "ZXDown", analysis, production, category)
        thetemplatesfiles += [ZXUp, ZXDown]

    if config.applyMINLOsystematics:
        MINLOUp = TemplatesFile(channel, "ggh", "MINLOUp", analysis, production, category)
        MINLODn = TemplatesFile(channel, "ggh", "MINLODn", analysis, production, category)
        thetemplatesfiles += [MINLOUp, MINLODn]

    outfiles = {templatesfile: ROOT.TFile(templatesfile.templatesfile(), "RECREATE") for templatesfile in thetemplatesfiles}

    store = []

    if config.applym4lshapesystematics:
        ggHuntaggedSM = {systematic:
                           gettemplate(channel, "ggH", analysis, production, "Untagged", analysis.purehypotheses[0], systematic)
                           .ProjectionZ().Clone("projection_{}".format(systematic))
                         for systematic in ("", "ScaleUp", "ScaleDown", "ResUp", "ResDown")}
        ggHnominal = ggHuntaggedSM[""]
        for syst in "ScaleUp", "ScaleDown", "ResUp", "ResDown":
            ggHsyst = ggHuntaggedSM[syst]
            for _ in "ggh", "vbf", "zh", "wh", "tth":
                if _ == "ggh" and category == "Untagged": continue
                tf = TemplatesFile(channel, _, analysis, production, category)
                for t in tf.templates() + tf.inttemplates():
                    h = t.gettemplate().Clone()
                    integral = h.Integral()
                    h.SetDirectory(outfiles[ScaleAndRes[_,syst]])
                    for x, y, z in izip(xrange(1, h.GetNbinsX()+1), xrange(1, h.GetNbinsY()+1), xrange(1, h.GetNbinsZ()+1)):
                        h.SetBinContent(
                                        x, y, z,
                                          h.GetBinContent(x, y, z)
                                          * ggHsyst.GetBinContent(z)
                                          / ggHnominal.GetBinContent(z)
                                       )
                    store.append(h)

        for templatesfile in ScaleAndRes.values(): outfiles[templatesfile].Write()

    if config.applym4lshapesystematics and config.combinem4lshapesystematics:
        for _ in "ggh", "vbf", "zh", "wh", "tth":
            tf = TemplatesFile(channel, _, analysis, production, category)
            for t in tf.templates() + tf.inttemplates():
                if isinstance(t, Template): hypothesis = t.hypothesis
                elif isinstance(t, IntTemplate): hypothesis = t.interferencetype
                h = gettemplate(t.productionmode, hypothesis, channel, analysis, production, category)
                #ScaleResDown = nominal + (ScaleDown-nominal) + (ResDown-nominal)
                #           = -nominal + ScaleDown + ResDown
                hdn = h.Clone()
                hdn.SetDirectory(outfiles[ScaleResDown[_]])
                hdn.Scale(-1)
                hdn.Add(gettemplate(t.productionmode, hypothesis, channel, analysis, "ResDown", production, category))
                hdn.Add(gettemplate(t.productionmode, hypothesis, channel, analysis, "ScaleDown", production, category))

                #ScaleResUp = nominal - (ScaleResDown-nominal)
                #             = 2*nominal - ScaleResDown
                hup = h.Clone()
                hup.SetDirectory(outfiles[ScaleResUp[_]])
                hup.Scale(2)
                hup.Add(hdn, -1)

                store += [hup, hdn]

    if config.applyMINLOsystematics and category != "Untagged":
        tf = TemplatesFile(channel, "ggh", analysis, production, category)
        SM = analysis.purehypotheses[0]
        POWHEG3DSM = gettemplate("ggH", SM, tf)
        POWHEG2DSM = POWHEG3DSM.Project3D("yxe")
        MINLO2DSM = POWHEG2DSM.Clone("MINLO2DSM")
        MINLO2DSM.Reset("M")
        for c in channels:
            MINLO2DSM.Add(gettemplate("ggH", "0+", "MINLO_SM", "MINLO", c, analysis, production, category))

        for t in tf.templates():
            POWHEG3D = t.gettemplate()
            MINLO3D = POWHEG3D.Clone(Template("ggH", t.hypothesis, "MINLOUp", channel, analysis, production, category).templatename())
            MINLO3D.SetDirectory(outfiles[MINLOUp])
            MINLO3DDn = POWHEG3D.Clone(Template("ggH", t.hypothesis, "MINLODn", channel, analysis, production, category).templatename())
            MINLO3DDn.SetDirectory(outfiles[MINLODn])
            for x, y, z in izip(xrange(1, MINLO3D.GetNbinsX()+1), xrange(1, MINLO3D.GetNbinsY()+1), xrange(1, MINLO3D.GetNbinsZ()+1)):
                MINLO3D.SetBinContent(
                                      x, y, z,
                                      MINLO3D.GetBinContent(x, y, z) * MINLO2DSM.GetBinContent(x, y) / POWHEG2DSM.GetBinContent(x, y)
                                     )
                MINLO3DDn.SetBinContent(
                                        x, y, z,
                                        POWHEG3D.GetBinContent(x, y, z) **2 / MINLO3D.GetBinContent(x, y, z)
                                       )
            store += [MINLO3D, MINLO3DDn]
        for templatesfile in MINLOUp, MINLODn: outfiles[templatesfile].Write()

        for t in tf.inttemplates():
            MINLOintUp = POWHEG3D.Clone(IntTemplate("ggH", t.interferencetype, "MINLOUp", channel, analysis, production, category).templatename())
            MINLOintUp.SetDirectory(outfiles[MINLOUp])
            MINLOintUp.Reset("M")
            for template, factor in t.templatesandfactors:
                MINLOintUp.Add(getattr(outfiles[MINLOUp], template.templatename()), factor)

            MINLOintDn = POWHEG3D.Clone(IntTemplate("ggH", t.interferencetype, "MINLODn", channel, analysis, production, category).templatename())
            MINLOintDn.SetDirectory(outfiles[MINLODn])
            MINLOintDn.Reset("M")
            for template, factor in t.templatesandfactors:
                MINLOintDn.Add(getattr(outfiles[MINLODn], template.templatename()), factor)

            store += [MINLOintUp, MINLOintDn]

    if config.applyZXshapesystematics:
        ZXtemplate = gettemplate("ZX", channel, analysis, production, category)
        qqZZtemplate = gettemplate("qqZZ", channel, analysis, production, category)
        ZXUptemplate = qqZZtemplate.Clone(Template("ZX", channel, analysis, "ZXUp", production, category).templatename())
        ZXUptemplate.Scale(ZXtemplate.Integral() / ZXUptemplate.Integral())
        ZXUptemplate.SetDirectory(outfiles[ZXUp])

        ZXDowntemplate = ZXtemplate.Clone(Template("ZX", channel, analysis, "ZXDown", production, category).templatename())
        ZXDowntemplate.SetDirectory(outfiles[ZXDown])
        ZXDowntemplate.Scale(2)
        ZXDowntemplate.Add(ZXUptemplate, -1)
        ZXDowntemplate.Floor()
        ZXDowntemplate.Scale(ZXtemplate.Integral() / ZXDowntemplate.Integral())

    for f in outfiles.values():
        f.Write()
        f.Close()

if __name__ == "__main__":
    for production in productions:
        for analysis in analyses:
            for channel in channels:
                for category in categories:
                    print production, analysis, channel, category
                    combinesystematics(channel, analysis, production, category)
