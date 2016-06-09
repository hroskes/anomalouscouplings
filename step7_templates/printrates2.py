from helperstuff.enums import Channel
import ROOT
luminosity = 25 #why not?

def printrates(flavor):
    flavor = Channel(flavor)
    f = ROOT.TFile(flavor.templatesfile(False))
    print f.template0PlusAdapSmoothMirror.Integral()*luminosity / 6,
    f = ROOT.TFile(flavor.templatesfile(True))
    print f.templateqqZZAdapSmoothMirror.Integral()*luminosity,
    print f.templateggZZAdapSmoothMirror.Integral()*luminosity,
    print f.templateZXAdapSmoothMirror.Integral()*luminosity,

if __name__ == "__main__":
    for flavor in "2e2mu", "4e", "4mu":
        print flavor
        print
        printrates(flavor)
        print
        print
