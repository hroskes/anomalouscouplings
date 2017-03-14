#!/usr/bin/env python

"""
combine ScaleUp, ResUp --> ScaleResUp
ZX systematics templates from qqZZ
"""

from helperstuff.enums import analyses, channels, productions, Template, TemplatesFile
from helperstuff.filemanager import tfiles
import ROOT

def combinesystematics(flavor, analysis, production):
    ScaleResUp = TemplatesFile(flavor, "signal", "ScaleResUp", analysis, production)
    ScaleResDown = TemplatesFile(flavor, "signal", "ScaleResDown", analysis, production)
    ZXUp = TemplatesFile(flavor, "bkg", "ZXUp", analysis, production)
    ZXDown = TemplatesFile(flavor, "bkg", "ZXDown", analysis, production)

    outfiles = {templatesfile: ROOT.TFile(templatesfile.templatesfile(), "RECREATE") for templatesfile in (ScaleResUp, ScaleResDown, ZXUp, ZXDown)}

    store = []

    for sample in analysis.signalsamples():
        h = Template(sample, flavor, analysis, production).gettemplate()
        #ScaleResUp = nominal + (ScaleUp-nominal) + (ResUp-nominal)
        #           = -nominal + ScaleUp + ResUp
        hup = h.Clone(Template(sample, flavor, analysis, "ScaleResUp", production).templatename())
        hup.SetDirectory(outfiles[ScaleResUp])
        hup.Scale(-1)
        hup.Add(Template(sample, flavor, analysis, "ResUp", production).gettemplate())
        hup.Add(Template(sample, flavor, analysis, "ScaleUp", production).gettemplate())

        #ScaleResDown = nominal - (ScaleResUp-nominal)
        #             = 2*nominal - ScaleResUp
        hdn = h.Clone(Template(sample, flavor, analysis, "ScaleResDown", production).templatename())
        hdn.SetDirectory(outfiles[ScaleResDown])
        hdn.Scale(2)
        hdn.Add(hup, -1)

        store += [hup, hdn]

    ZXtemplate = Template("ZX", flavor, analysis, production).gettemplate()
    qqZZtemplate = Template("qqZZ", flavor, analysis, production).gettemplate()
    ZXUptemplate = qqZZtemplate.Clone(Template("ZX", flavor, analysis, "ZXUp", production).templatename())
    ZXUptemplate.Scale(ZXtemplate.Integral() / ZXUptemplate.Integral())
    ZXUptemplate.SetDirectory(outfiles[ZXUp])

    ZXDowntemplate = ZXtemplate.Clone(Template("ZX", flavor, analysis, "ZXDown", production).templatename())
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
                print production, analysis, channel
                combinesystematics(channel, analysis, production)
