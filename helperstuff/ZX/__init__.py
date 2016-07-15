import os
import ROOT
ZXfolder = os.path.dirname(__file__)
ReducibleBackgroundFile = os.path.join(ZXfolder, "ReducibleBackgroundAA_2015.C")
ROOT.gROOT.LoadMacro(ReducibleBackgroundFile+"+")

def setup(production):
    release = int(production.release)
    assert release in (76, 80)
    ROOT.setup(release, ZXfolder)
