import os
import ROOT
cconstantsfolder = os.path.dirname(__file__)
cconstantsfile = os.path.join(cconstantsfolder, "cconstants.C")
ROOT.gROOT.LoadMacro(cconstantsfile+"+")

getDVBF2jetsConstant = ROOT.getDVBF2jetsConstant
getDVBF1jetConstant = ROOT.getDVBF1jetConstant
getDbkgkinConstant = ROOT.getDbkgkinConstant
getDbkgConstant = ROOT.getDbkgConstant
