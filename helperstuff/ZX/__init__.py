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
#https://github.com/CJLST/ZZAnalysis/blob/b4cc949af2ca81a9dd2bfeaba47768c7ea0cfd13/AnalysisStep/data/FakeRates/FakeRate_SS_2016D.root
convertTGraphstoTH1Fs.convertTGraphstoTH1Fs(os.path.join(ZXfolder, "FakeRate_SS_2016D.root"))
#https://github.com/CJLST/ZZAnalysis/blob/4b0ac1ccda2a60295e8233069dcd1a17802894f7/AnalysisStep/data/FakeRates/FakeRate_SS_2016D.root
convertTGraphstoTH1Fs.convertTGraphstoTH1Fs(os.path.join(ZXfolder, "FakeRate_SS_2016D_12.9fb-1.root"))

from ROOT import CRZLLss, fakeRate13TeV, test_bit
