import os
import ROOT
folder = os.path.dirname(__file__)
file = os.path.join(folder, "TextExtentTest.C")
ROOT.gROOT.LoadMacro(file+"+")

GetTextWidth = ROOT.GetTextWidth
