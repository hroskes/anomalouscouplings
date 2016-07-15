import os
import ROOT
ZXfolder = os.path.dirname(__file__)
ReducibleBackgroundFile = os.path.join(folder, "ReducibleBackgroundAA_2015.C")
ROOT.gROOT.LoadMacro(ReducibleBackgroundFile+"+")

def setup(release):
    release = int(release)
    assert release in (76, 80)
    ROOT.setup(release, folder)
