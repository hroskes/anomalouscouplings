"""
combine ScaleUp, ResUp --> ScaleResUp
ZX systematics templates from qqZZ
"""

from helperstuff.enums import analyses, channels, Template, TemplatesFile
from helperstuff.filemanager import tfiles
import ROOT

def combinesystematics(flavor, analysis):
    ScaleResUp = TemplatesFile(flavor, "signal", "ScaleResUp", analysis)
    ScaleResDown = TemplatesFile(flavor, "signal", "ScaleResDown", analysis)
    ZXUp = TemplatesFile(flavor, "bkg", "ZXUp", analysis)
    ZXDown = TemplatesFile(flavor, "bkg", "ZXDown", analysis)

    outfiles = {templatesfile: ROOT.TFile(templatesfile.templatesfile(), "RECREATE") for templatesfile in (ScaleResUp, ScaleResDown, ZXUp, ZXDown)}

    store = []

    for sample in analysis.signalsamples():
        h = Template(sample, flavor, analysis).gettemplate()
        #ScaleResUp = nominal + (ScaleUp-nominal) + (ResUp-nominal)
        #           = -nominal + ScaleUp + ResUp
        hup = h.Clone(Template(sample, flavor, analysis, "ScaleResUp").templatename())
        hup.SetDirectory(outfiles[ScaleResUp])
        hup.Scale(-1)
        hup.Add(Template(sample, flavor, analysis, "ResUp").gettemplate())
        hup.Add(Template(sample, flavor, analysis, "ScaleUp").gettemplate())

        #ScaleResDown = nominal - (ScaleResUp-nominal)
        #             = 2*nominal - ScaleResUp
        hdn = h.Clone(Template(sample, flavor, analysis, "ScaleResDown").templatename())
        hdn.SetDirectory(outfiles[ScaleResDown])
        hdn.Scale(2)
        hdn.Add(hup, -1)

        store += [hup, hdn]

    ZXtemplate = Template("ZX", flavor, analysis).gettemplate()
    qqZZtemplate = Template("qqZZ", flavor, analysis).gettemplate()
    ZXUptemplate = qqZZtemplate.Clone(Template("ZX", flavor, analysis, "ZXUp").templatename())
    ZXUptemplate.Scale(ZXtemplate.Integral() / ZXUptemplate.Integral())
    ZXUptemplate.SetDirectory(outfiles[ZXUp])

    ZXDowntemplate = ZXtemplate.Clone(Template("ZX", flavor, analysis, "ZXDown").templatename())
    ZXDowntemplate.SetDirectory(outfiles[ZXDown])
    ZXDowntemplate.Scale(2)
    ZXDowntemplate.Add(ZXUptemplate, -1)

    print flavor, ZXtemplate.Integral(), ZXUptemplate.Integral(), ZXDowntemplate.Integral()

    for f in outfiles.values():
        f.Write()
        f.Close()

if __name__ == "__main__":
    for analysis in analyses:
        for channel in channels:
            combinesystematics(channel, analysis)
