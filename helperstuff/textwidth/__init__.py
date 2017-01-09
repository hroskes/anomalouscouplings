import os
import ROOT
import utilities
folder = os.path.dirname(__file__)
file = os.path.join(folder, "TextExtentTest.C")
utilities.LoadMacro(file+"+")

GetTextWidth = ROOT.GetTextWidth
