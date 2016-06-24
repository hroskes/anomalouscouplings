from helperstuff.enums import channels, TemplatesFile
from helperstuff.filemanager import tfiles
import ROOT

"""combine ScaleUp, ResUp --> ScaleResUp"""

def combinesystematics(flavor):
    nominal = TemplatesFile(flavor, "signal")
    ResUp = TemplatesFile(flavor, "signal", "ResUp")
    ScaleUp = TemplatesFile(flavor, "signal", "ScaleUp")
    ScaleResUp = TemplatesFile(flavor, "signal", "ScaleResUp")
    ScaleResDown = TemplatesFile(flavor, "signal", "ScaleResDown")

    infiles = {templatesfile: tfiles[templatesfile.templatesfile()] for templatesfile in (nominal, ScaleUp, ResUp)}
    outfiles = {templatesfile: ROOT.TFile(templatesfile.templatesfile(), "RECREATE") for templatesfile in (ScaleResUp, ScaleResDown)}

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

        print flavor, h.GetName(), infiles[ResUp].Get(name).Integral() + infiles[ScaleUp].Get(name).Integral() - h.Integral() - hup.Integral()

    for f in outfiles.values():
        f.Write()
        f.Close()

if __name__ == "__main__":
    for channel in channels:
        combinesystematics(channel)
