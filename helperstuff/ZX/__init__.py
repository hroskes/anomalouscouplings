import os
import ROOT
ZXfolder = os.path.dirname(__file__)
ReducibleBackgroundFile = os.path.join(ZXfolder, "ReducibleBackgroundAA_2015.C")
ROOT.gROOT.LoadMacro(ReducibleBackgroundFile+"+")

def setup(production):
    success = ROOT.setup(int(production), ZXfolder)
    if not success:
        raise ValueError("Bad production: {}".format(production))

import convertTGraphstoTH1Fs
convertTGraphstoTH1Fs.convertTGraphstoTH1Fs(os.path.join(ZXfolder, "FakeRate_SS_2016B.root"))
convertTGraphstoTH1Fs.convertTGraphstoTH1Fs(os.path.join(ZXfolder, "FakeRate_SS_2016D.root"))

from ROOT import CRZLLss, fakeRate13TeV, test_bit
