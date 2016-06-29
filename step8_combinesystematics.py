"""
combine ScaleUp, ResUp --> ScaleResUp
ZX systematics templates from qqZZ
"""

from helperstuff.enums import analyses, channels, TemplatesFile
from helperstuff.filemanager import tfiles
import ROOT

def combinesystematics(flavor, analysis):
    nominal = TemplatesFile(flavor, "signal", analysis)
    ResUp = TemplatesFile(flavor, "signal", "ResUp", analysis)
    ScaleUp = TemplatesFile(flavor, "signal", "ScaleUp", analysis)
    ScaleResUp = TemplatesFile(flavor, "signal", "ScaleResUp", analysis)
    ScaleResDown = TemplatesFile(flavor, "signal", "ScaleResDown", analysis)
    bkg = TemplatesFile(flavor, "bkg", analysis)
    ZXUp = TemplatesFile(flavor, "bkg", "ZXUp", analysis)
    ZXDown = TemplatesFile(flavor, "bkg", "ZXDown", analysis)

    infiles = {templatesfile: tfiles[templatesfile.templatesfile()] for templatesfile in (nominal, ScaleUp, ResUp, bkg)}
    outfiles = {templatesfile: ROOT.TFile(templatesfile.templatesfile(), "RECREATE") for templatesfile in (ScaleResUp, ScaleResDown, ZXUp, ZXDown)}

    objects = [key.ReadObj() for key in infiles[nominal].GetListOfKeys()]
    hists = [h for h in objects if h.GetName() != "controlPlots"]

    store = []

    for h in hists:
        if not isinstance(h, ROOT.TH1):
            raise TypeError("Unknown thing {} of type {} in root file!".format(h, type(h)))
        #ScaleResUp = nominal + (ScaleUp-nominal) + (ResUp-nominal)
        #           = -nominal + ScaleUp + ResUp
        hup = h.Clone()
        hup.SetDirectory(outfiles[ScaleResUp])
        hup.Scale(-1)
        name = hup.GetName()
        hup.Add(infiles[ResUp].Get(name))
        hup.Add(infiles[ScaleUp].Get(name))

        #ScaleResDown = nominal - (ScaleResUp-nominal)
        #             = 2*nominal - ScaleResUp
        hdn = h.Clone()
        hdn.SetDirectory(outfiles[ScaleResDown])
        hdn.Scale(2)
        hdn.Add(hup, -1)

        store += [hup, hdn]

    if analysis.domirror():
        ZXtemplate = infiles[bkg].templateZXAdapSmoothMirror
        qqZZtemplate = infiles[bkg].templateqqZZAdapSmoothMirror
    else:
        ZXtemplate = infiles[bkg].templateZXAdapSmooth
        qqZZtemplate = infiles[bkg].templateqqZZAdapSmooth
    ZXUptemplate = qqZZtemplate.Clone(ZXtemplate.GetName())
    ZXUptemplate.Scale(ZXtemplate.Integral() / ZXUptemplate.Integral())
    ZXUptemplate.SetDirectory(outfiles[ZXUp])

    ZXDowntemplate = ZXtemplate.Clone()
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
