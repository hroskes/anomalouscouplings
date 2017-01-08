import utilities
import downloadfromCJLST
import os
import ROOT

CJLSTscriptsfolder = os.path.dirname(os.path.abspath(__file__))
CJLSTscriptsfile = os.path.join(CJLSTscriptsfolder, "CJLSTscripts.cc")

#have to be in order of who includes whose header file
scripts = ["cConstants", "Category", "bitops", "FinalStates"]

downloader = downloadfromCJLST.Downloader("01316520484cc2c15443ab6f3419d59b20b1a715")
for script in scripts:
    downloader.add("AnalysisStep/src/{}.cc".format(script))
    downloader.add("AnalysisStep/interface/{}.h".format(script))

downloader.add("AnalysisStep/test/Macros/ReducibleBackgroundAA_2015.C")
for rootfile in "FakeRate_SS_2016D.root",:
    downloader.add(os.path.join("AnalysisStep/data/FakeRates", rootfile))

with utilities.cd(CJLSTscriptsfolder):
    downloader.download()

for script in scripts:
    ROOT.gROOT.LoadMacro(os.path.join(CJLSTscriptsfolder, script+".cc+"))
ROOT.gROOT.LoadMacro(os.path.join(CJLSTscriptsfolder, "ReducibleBackgroundAA_2015_adapted.C+"))

from ROOT import categoryIchep16, getDVBF2jetsConstant, getDVBF1jetConstant, getDbkgkinConstant, getDbkgConstant, UntaggedIchep16, VBF1jTaggedIchep16, VBF2jTaggedIchep16, VHLeptTaggedIchep16, VHHadrTaggedIchep16, ttHTaggedIchep16
