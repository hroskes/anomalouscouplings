import utilities
import downloadfromCJLST
import os
import ROOT

CJLSTscriptsfolder = os.path.dirname(os.path.abspath(__file__))

#have to be in order of who includes whose header file
scripts = ["cConstants", "Discriminants", "Category", "bitops", "FinalStates"]

downloader = downloadfromCJLST.Downloader("19cadfb4e01a68d9f14c6697c1ca7f3d9d339bbc")
for script in scripts:
    downloader.add("AnalysisStep/src/{}.cc".format(script))
    downloader.add("AnalysisStep/interface/{}.h".format(script))

downloader.add("AnalysisStep/test/Macros/ReducibleBackgroundAA_2015.C")
for rootfile in "FakeRate_SS_Moriond368.root",:
    downloader.add(os.path.join("AnalysisStep/data/FakeRates", rootfile))
for cconstant in "Dbkgkin_2e2mu", "Dbkgkin_4e", "Dbkgkin_4mu", "DjVBF", "DjjVBF", "DjjWH", "DjjZH":
    downloader.add("AnalysisStep/data/cconstants/SmoothKDConstant_m4l_{}13TeV.root".format(cconstant))

with utilities.cd(CJLSTscriptsfolder):
    downloader.download()

for script in scripts:
    utilities.LoadMacro(os.path.join(CJLSTscriptsfolder, script+".cc+"))
utilities.LoadMacro(os.path.join(CJLSTscriptsfolder, "ReducibleBackgroundAA_2015_adapted.C+"))

from ROOT import categoryMor17, UntaggedMor17, VBF1jTaggedMor17, VBF2jTaggedMor17, VHLeptTaggedMor17, VHHadrTaggedMor17, ttHTaggedMor17, VHMETTaggedMor17

from ROOT import getDVBF2jetsConstant, getDVBF1jetConstant, getDWHhConstant, getDZHhConstant, getDbkgkinConstant, getDbkgConstant
from ROOT import getDVBF2jetsWP, getDVBF1jetWP, getDWHhWP, getDZHhWP
from ROOT import getDVBF2jetsConstant_shiftWP, getDVBF1jetConstant_shiftWP, getDWHhConstant_shiftWP, getDZHhConstant_shiftWP
