import utilities
import downloadfromCJLST
import os
import ROOT

CJLSTscriptsfolder = os.path.dirname(__file__)
CJLSTscriptsfile = os.path.join(CJLSTscriptsfolder, "CJLSTscripts.cc")

#have to be in order of who includes whose header file
scripts = ["cConstants", "Category"]

with utilities.cd(CJLSTscriptsfolder):
    for script in scripts:
        downloadfromCJLST.download("AnalysisStep/src/{}.cc".format(script))
        downloadfromCJLST.download("AnalysisStep/interface/{}.h".format(script))

for script in scripts:
    ROOT.gROOT.LoadMacro(os.path.join(CJLSTscriptsfolder, script+".cc+"))

from ROOT import categoryIchep16, getDVBF2jetsConstant, getDVBF1jetConstant, getDbkgkinConstant, getDbkgConstant, UntaggedIchep16, VBF1jTaggedIchep16, VBF2jTaggedIchep16, VHLeptTaggedIchep16, VHHadrTaggedIchep16, ttHTaggedIchep16
