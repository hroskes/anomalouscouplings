from helperstuff.enums import Channel
import ROOT
luminosity = 1 #why not?

def printrates(flavor):
    flavor = Channel(flavor)
    f = ROOT.TFile(flavor.templatesfile(False))
    print "rate ggH", f.template0PlusAdapSmoothMirror.Integral()*luminosity
    f = ROOT.TFile(flavor.templatesfile(True))
    print "rate qqZZ", f.templateqqZZAdapSmoothMirror.Integral()*luminosity
    print "rate ggZZ", f.templateggZZAdapSmoothMirror.Integral()*luminosity
    print "rate zjets", f.templateZXAdapSmoothMirror.Integral()*luminosity

if __name__ == "__main__":
    for flavor in "2e2mu", "4e", "4mu":
        print flavor
        print
        printrates(flavor)
        print
