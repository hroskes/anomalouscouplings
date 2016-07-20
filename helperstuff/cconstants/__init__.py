import os
import ROOT
cconstantsfolder = os.path.dirname(__file__)
cconstantsfile = os.path.join(cconstantsfolder, "cconstants.C")
ROOT.gROOT.LoadMacro(cconstantsfile+"+")

from ROOT import getDVBF2jetsConstant, getDVBF1jetConstant, getDbkgkinConstant, getDbkgConstant
