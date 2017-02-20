#!/usr/bin/env python
"""
combine ScaleUp, ResUp --> ScaleResUp
ZX systematics templates from qqZZ
"""

from helperstuff import config
from helperstuff.enums import analyses, categories, channels, productions
from helperstuff.templates import Template, TemplatesFile
from helperstuff.utilities import tfiles
import ROOT

def combinesystematics(channel, analysis, production, category):
    thetemplatesfiles = []

    if config.applym4lshapesystematics and config.combinem4lshapesystematics:
        ScaleResUp = TemplatesFile(channel, "signal", "ScaleResUp", analysis, production, category)
        ScaleResDown = TemplatesFile(channel, "signal", "ScaleResDown", analysis, production, category)
        thetemplatesfiles += [ScaleResUp, ScaleResDown]

    if config.applyZXshapesystematics:
        ZXUp = TemplatesFile(channel, "bkg", "ZXUp", analysis, production, category)
        ZXDown = TemplatesFile(channel, "bkg", "ZXDown", analysis, production, category)
        thetemplatesfiles += [ZXUp, ZXDown]

    outfiles = {templatesfile: ROOT.TFile(templatesfile.templatesfile(), "RECREATE") for templatesfile in thetemplatesfiles}

    store = []

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
