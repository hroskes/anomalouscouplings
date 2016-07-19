from array import array
from math import sqrt
import ROOT

def convertTGraphtoTH1F(g):
    n = g.GetN()
    if isinstance(g, ROOT.TGraphAsymmErrors):
        X, Y, eXup, eXdn, eYup, eYdn = [[z for i, z in zip(range(n), thing)] for thing in g.GetX(), g.GetY(), g.GetEXhigh(), g.GetEXlow(), g.GetEYhigh(), g.GetEYlow()]
    elif isinstance(g, ROOT.TGraphErrors):
        X, Y, eXup, eXdn, eYup, eYdn = [[z for i, z in zip(range(n), thing)] for thing in g.GetX(), g.GetY(), g.GetEX(), g.GetEX(), g.GetEY(), g.GetEY()]
    else:
        assert False
    xarray = [X[0]]
    for x, exup in zip(X, eXup):
        xarray.append(x+exup)
    print xarray
    h = ROOT.TH1F(g.GetName(), g.GetTitle(), n, array("d", xarray))
    for i, y, eyup, eydn in zip(range(1, n+1), Y, eYup, eYdn):
        h.SetBinContent(i, y)
        h.SetBinError(i, sqrt((eyup**2 + eydn**2)/2))
    return h

def convertTGraphstoTH1Fs(filename)
    if not os.path.exists(filename): raise OSError("{} does not exist!".format(filename))
    newfilename = f.replace(".root", "_hists.root")
    if os.path.exists(newfilename): return
    f = ROOT.TFile(filename)
    newf = ROOT.TFile(newfilename, "RECREATE")
    hists = [convertTGraphtoTH1F(k.ReadObj()) for k in f.GetListOfKeys()]
    newf.Write()
