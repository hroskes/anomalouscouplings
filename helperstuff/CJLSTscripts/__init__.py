import utilities
import downloadfromCJLST
import os
import ROOT

CJLSTscriptsfolder = os.path.dirname(__file__)
if CJLSTscriptsfolder == "": CJLSTscriptsfolder = "."
CJLSTscriptsfile = os.path.join(CJLSTscriptsfolder, "CJLSTscripts.cc")

#have to be in order of who includes whose header file
scripts = ["cConstants", "Category"]

with utilities.cd(CJLSTscriptsfolder):
    downloader = downloadfromCJLST.Downloader("01316520484cc2c15443ab6f3419d59b20b1a715")
    for script in scripts:
        downloader.add("AnalysisStep/src/{}.cc".format(script))
        downloader.add("AnalysisStep/interface/{}.h".format(script))
    downloader.download()

for script in scripts:
    ROOT.gROOT.LoadMacro(os.path.join(CJLSTscriptsfolder, script+".cc+"))

from ROOT import categoryIchep16, getDVBF2jetsConstant, getDVBF1jetConstant, getDbkgkinConstant, getDbkgConstant, UntaggedIchep16, VBF1jTaggedIchep16, VBF2jTaggedIchep16, VHLeptTaggedIchep16, VHHadrTaggedIchep16, ttHTaggedIchep16
