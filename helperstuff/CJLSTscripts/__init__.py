import utilities
import downloadfromCJLST
import os
import ROOT

CJLSTscriptsfolder = os.path.dirname(os.path.abspath(__file__))

#have to be in order of who includes whose header file
scripts = ["cConstants", "Discriminants", "Category", "bitops", "FinalStates"]

downloader = downloadfromCJLST.Downloader("15ba428aef070653bb720673b2f250c589258565")
for script in scripts:
    downloader.add("AnalysisStep/src/{}.cc".format(script))
    downloader.add("AnalysisStep/interface/{}.h".format(script))

downloader.add("AnalysisStep/test/ZpXEstimation/include/FakeRates.h")
downloader.add("AnalysisStep/test/ZpXEstimation/src/FakeRates.cpp")
for cconstant in "Dbkgkin_2e2mu", "Dbkgkin_4e", "Dbkgkin_4mu", "DjVBF", "DjjVBF", "DjjWH", "DjjZH", "DbkgjjEWQCD_4l_HadVHTagged_", "DbkgjjEWQCD_4l_JJVBFTagged_", "DbkgjjEWQCD_2l2l_HadVHTagged_", "DbkgjjEWQCD_2l2l_JJVBFTagged_":
    downloader.add("AnalysisStep/data/cconstants/SmoothKDConstant_m4l_{}13TeV.root".format(cconstant))
for rootfile in "FakeRates_SS_2016_Legacy.root", "FakeRates_SS_2017_Legacy.root", "FakeRates_SS_2018_Legacy.root":
    downloader.add(os.path.join("AnalysisStep/data/FakeRates", rootfile))
for rootfile in "FakeRate_SS_Moriond368.root", "FakeRates_SS_Moriond18.root", "FakeRates_SS_Moriond19.root":
    downloader.add(os.path.join("AnalysisStep/data/FakeRates", rootfile), sha1="e280b6b44d768d52602b7475f0bf724d3aca8531^")

with utilities.cd(CJLSTscriptsfolder):
    downloader.download()

for script in scripts:
    utilities.LoadMacro(os.path.join(CJLSTscriptsfolder, script+".cc+"))
utilities.LoadMacro(os.path.join(CJLSTscriptsfolder, "FakeRates.cpp+"))

from ROOT import categoryAC19, UntaggedAC19, VBF1jTaggedAC19, VBF2jTaggedAC19, VHLeptTaggedAC19, VHHadrTaggedAC19, ttHLeptTaggedAC19, ttHHadrTaggedAC19, VHMETTaggedAC19, BoostedAC19, categoryMor18

from ROOT import getDVBF2jetsConstant, getDVBF1jetConstant, getDWHhConstant, getDZHhConstant, getDbkgkinConstant, getDbkgConstant
from ROOT import getDVBF2jetsWP, getDVBF1jetWP, getDWHhWP, getDZHhWP
from ROOT import getDVBF2jetsConstant_shiftWP, getDVBF1jetConstant_shiftWP, getDWHhConstant_shiftWP, getDZHhConstant_shiftWP

from ROOT import D_bkg_VBFdec, D_bkg_VHdec, DVBF1j_ME
