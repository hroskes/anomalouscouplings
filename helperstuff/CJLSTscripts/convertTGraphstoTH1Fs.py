from array import array
from math import sqrt
import os

import ROOT

from utilities import OneAtATime

def convertTGraphtoTH1F(key):
    g = key.ReadObj()
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
    h = ROOT.TH1F(key.GetName(), key.GetTitle(), n, array("d", xarray))
    print h.GetName(), xarray
    for i, y, eyup, eydn in zip(range(1, n+1), Y, eYup, eYdn):
        h.SetBinContent(i, y)
        h.SetBinError(i, sqrt((eyup**2 + eydn**2)/2))
    return h

def convertTGraphstoTH1Fs(filename):
    if not os.path.exists(filename): raise OSError("{} does not exist!".format(filename))
    newfilename = filename.replace(".root", "_hists.root")
    with OneAtATime(newfilename+".tmp", 2, task="converting TGraphs to TH1Fs"):
        if os.path.exists(newfilename): return
        f = ROOT.TFile(filename)
        try:
            newf = ROOT.TFile(newfilename, "RECREATE")
            hists = [convertTGraphtoTH1F(k) for k in f.GetListOfKeys()]
            newf.Write()
        except:
            try:
                os.remove(newfilename)
            except:
                pass
            raise
