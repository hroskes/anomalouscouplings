#!/usr/bin/env python
"""
combine ScaleUp, ResUp --> ScaleResUp
ZX systematics templates from qqZZ
"""
from itertools import izip

import ROOT

from helperstuff import config
from helperstuff.combinehelpers import gettemplate
from helperstuff.enums import analyses, categories, channels, productions
from helperstuff.templates import Template, TemplatesFile
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

    if config.applym4lshapesystematics and config.combinem4lshapesystematics:
        for sample in analysis.signalsamples():
            h = Template(sample, channel, analysis, production, category).gettemplate()
            #ScaleResUp = nominal + (ScaleUp-nominal) + (ResUp-nominal)
            #           = -nominal + ScaleUp + ResUp
            hup = h.Clone(Template(sample, channel, analysis, "ScaleResUp", production, category).templatename())
            hup.SetDirectory(outfiles[ScaleResUp])
            hup.Scale(-1)
            hup.Add(Template(sample, channel, analysis, "ResUp", production, category).gettemplate())
            hup.Add(Template(sample, channel, analysis, "ScaleUp", production, category).gettemplate())

            #ScaleResDown = nominal - (ScaleResUp-nominal)
            #             = 2*nominal - ScaleResUp
            hdn = h.Clone(Template(sample, channel, analysis, "ScaleResDown", production, category).templatename())
            hdn.SetDirectory(outfiles[ScaleResDown])
            hdn.Scale(2)
            hdn.Add(hup, -1)

            store += [hup, hdn]

    if config.applyZXshapesystematics:
        ZXtemplate = Template("ZX", channel, analysis, production, category).gettemplate()
        qqZZtemplate = Template("qqZZ", channel, analysis, production, category).gettemplate()
        ZXUptemplate = qqZZtemplate.Clone(Template("ZX", channel, analysis, "ZXUp", production, category).templatename())
        ZXUptemplate.Scale(ZXtemplate.Integral() / ZXUptemplate.Integral())
        ZXUptemplate.SetDirectory(outfiles[ZXUp])

        ZXDowntemplate = ZXtemplate.Clone(Template("ZX", channel, analysis, "ZXDown", production, category).templatename())
        ZXDowntemplate.SetDirectory(outfiles[ZXDown])
        ZXDowntemplate.Scale(2)
        ZXDowntemplate.Add(ZXUptemplate, -1)

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
