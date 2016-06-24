import os
import ROOT
ReducibleBackgroundFile = os.path.join(os.path.dirname(__file__), "ReducibleBackgroundAA_2015.C")
ROOT.gROOT.LoadMacro(ReducibleBackgroundFile+"+")
