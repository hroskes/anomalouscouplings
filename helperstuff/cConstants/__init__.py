from .. import filemanager
import downloadfromCJLST
import os
import ROOT

cConstantsfolder = os.path.dirname(__file__)
cConstantsfile = os.path.join(cConstantsfolder, "cConstants.cc")

with filemanager.cd(cConstantsfolder):
    downloadfromCJLST.download("AnalysisStep/src/cConstants.cc")
    downloadfromCJLST.download("AnalysisStep/interface/cConstants.h")


ROOT.gROOT.LoadMacro(cConstantsfile+"+")

from ROOT import getDVBF2jetsConstant, getDVBF1jetConstant, getDbkgkinConstant, getDbkgConstant
