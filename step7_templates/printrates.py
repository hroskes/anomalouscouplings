import collections
from helperstuff.enums import Channel, channels
from helperstuff.filemanager import tfiles
import ROOT
import sys

luminosity = 10

def printrates(flavor, mode):
    flavor = Channel(flavor)
    f = tfiles[flavor.templatesfile(False)]
    ggH = f.template0PlusAdapSmoothMirror.Integral()*luminosity
    f = tfiles[flavor.templatesfile(True)]
    qqZZ = f.templateqqZZAdapSmoothMirror.Integral()*luminosity
    ggZZ = f.templateggZZAdapSmoothMirror.Integral()*luminosity
    ZX = f.templateZXAdapSmoothMirror.Integral() * (
         #no idea about the absolute
         #rescale to Moriond
             (0.408547+0.311745+0.0106453+0.716686+0.0199815)  #email from Simon, "Inputs for the cards", Feb 9 at 4:56AM
              / sum(tfiles[c.templatesfile(True)].templateZXAdapSmoothMirror.Integral() for c in channels)
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
        printrates(flavor, int(sys.argv[1]))
        print
