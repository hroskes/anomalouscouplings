import utilities
import downloadfromCJLST
import os
import ROOT

CJLSTscriptsfolder = os.path.dirname(os.path.abspath(__file__))

#have to be in order of who includes whose header file
scripts = ["cConstants", "Category", "bitops", "FinalStates"]

downloader = downloadfromCJLST.Downloader("16a57b39523211dca423b18d3fa968af9752b8ec")
for script in scripts:
    downloader.add("AnalysisStep/src/{}.cc".format(script))
    downloader.add("AnalysisStep/interface/{}.h".format(script))

downloader.add("AnalysisStep/test/Macros/ReducibleBackgroundAA_2015.C")
for rootfile in "FakeRate_SS_2016D.root",:
    downloader.add(os.path.join("AnalysisStep/data/FakeRates", rootfile))

with utilities.cd(CJLSTscriptsfolder):
    downloader.download()

for script in scripts:
    utilities.LoadMacro(os.path.join(CJLSTscriptsfolder, script+".cc+"))
utilities.LoadMacro(os.path.join(CJLSTscriptsfolder, "ReducibleBackgroundAA_2015_adapted.C+"))

from ROOT import categoryIchep16, UntaggedIchep16, VBF1jTaggedIchep16, VBF2jTaggedIchep16, VHLeptTaggedIchep16, VHHadrTaggedIchep16, ttHTaggedIchep16
from ROOT import getDVBF2jetsConstant, getDVBF1jetConstant, getDWHhConstant, getDZHhConstant, getDbkgkinConstant, getDbkgConstant
from ROOT import getDVBF2jetsWP, getDVBF1jetWP, getDWHhWP, getDZHhWP
from ROOT import getDVBF2jetsConstant_shiftWP, getDVBF1jetConstant_shiftWP, getDWHhConstant_shiftWP, getDZHhConstant_shiftWP
