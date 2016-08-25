from .. import filemanager
import downloadfromCJLST
import os
import ROOT

from .. import cConstants

Categoryfolder = os.path.dirname(__file__)
Categoryfile = os.path.join(Categoryfolder, "Category.cc")

with filemanager.cd(Categoryfolder):
    downloadfromCJLST.download("AnalysisStep/src/Category.cc")
    downloadfromCJLST.download("AnalysisStep/interface/Category.h")


ROOT.gROOT.LoadMacro(Categoryfile+"+")

from ROOT import categoryLegacy
