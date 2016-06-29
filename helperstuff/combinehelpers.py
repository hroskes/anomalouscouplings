import collections
from helperstuff import config
from helperstuff.enums import Channel, channels, TemplatesFile
from helperstuff.filemanager import tfiles
from helperstuff.samples import Sample
import ROOT

def getrates(flavor, analysis):
    flavor = Channel(flavor)
    f = tfiles[TemplatesFile(flavor, "signal", analysis).templatesfile()]
    ggH = f.template0PlusAdapSmoothMirror.Integral()*config.luminosity
    #other signal samples
    for productionmode in "VBF", "WplusH", "WminusH", "ZH", "ttH":
        sample = Sample(productionmode, "0+")
        f = tfiles[sample.withdiscriminantsfile()]
        t = f.candTree
        ZZFlav = flavor.ZZFlav()
        additionalxsec = 0
        for event in t:
            if 105 < t.ZZMass < 140 and t.Z1Flav*t.Z2Flav == ZZFlav:
                additionalxsec += getattr(t, sample.weightname())
        ggH += additionalxsec * config.luminosity


    f = tfiles[TemplatesFile(flavor, "bkg", analysis).templatesfile()]
    qqZZ = f.templateqqZZAdapSmoothMirror.Integral()*config.luminosity
    ggZZ = f.templateggZZAdapSmoothMirror.Integral()*config.luminosity
    ZX = f.templateZXAdapSmoothMirror.Integral() * (
         #no idea about the absolute
         #rescale to Moriond
             (0.408547+0.311745+0.0106453+0.716686+0.0199815)  #email from Simon, "Inputs for the cards", Feb 9 at 4:56AM
              / sum(tfiles[TemplatesFile(c, "bkg", analysis).templatesfile()].templateZXAdapSmoothMirror.Integral() for c in channels)
         # * ratio of luminosities
             * config.luminosity / 2.8
         )

    result =  "rate {} {} {} {}".format(ggH, qqZZ, ggZZ, ZX)
    return result
