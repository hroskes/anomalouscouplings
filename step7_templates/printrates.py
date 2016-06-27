import collections
from helperstuff.enums import Channel, channels, TemplatesFile
from helperstuff.filemanager import tfiles
from helperstuff.samples import Sample
import ROOT

luminosity = 10
mode = 2

def printrates(flavor):
    flavor = Channel(flavor)
    f = tfiles[TemplatesFile(flavor, "signal").templatesfile()]
    ggH = f.template0PlusAdapSmoothMirror.Integral()*luminosity
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
        ggH += additionalxsec * luminosity


    f = tfiles[TemplatesFile(flavor, "bkg").templatesfile()]
    qqZZ = f.templateqqZZAdapSmoothMirror.Integral()*luminosity
    ggZZ = f.templateggZZAdapSmoothMirror.Integral()*luminosity
    ZX = f.templateZXAdapSmoothMirror.Integral() * (
         #no idea about the absolute
         #rescale to Moriond
             (0.408547+0.311745+0.0106453+0.716686+0.0199815)  #email from Simon, "Inputs for the cards", Feb 9 at 4:56AM
              / sum(tfiles[TemplatesFile(c, "bkg").templatesfile()].templateZXAdapSmoothMirror.Integral() for c in channels)
         # * ratio of luminosities
             * luminosity / 2.8
         )

    if mode == 1:
        print "rate ggH", ggH
        print "rate qqZZ", qqZZ
        print "rate ggZZ", ggZZ
        print "rate zjets", ZX
    elif mode == 2:
        print "rate", ggH, qqZZ, ggZZ, ZX
    else:
        assert False

if __name__ == "__main__":
    for flavor in "2e2mu", "4e", "4mu":
        print flavor
        print
        printrates(flavor)
        print
