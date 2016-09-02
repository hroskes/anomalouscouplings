from helperstuff import config, style  #style needed to remove the statbox and title and adjust the axis label size
from helperstuff.enums import analyses, channels, productions, treesystematics
from helperstuff.templates import TemplatesFile
import os
import ROOT
import subprocess

def makecontrolplots(*args):
    templatesfile = TemplatesFile(*args)
    f = ROOT.TFile.Open(templatesfile.templatesfile())
    d = f.controlPlots

    split = os.path.split(templatesfile.templatesfile())
    saveasdir = os.path.join(config.plotsbasedir, "templateprojections", "controlplots", split[1].replace(".root", ""))
    try:
        os.makedirs(saveasdir)
    except OSError:
        pass

    exts = "png", "eps", "root", "pdf"

    axistitles = [templatesfile.analysis.purediscriminant(True), templatesfile.analysis.mixdiscriminant(True), templatesfile.systematic.D_bkg_0plus(True)]

    for key in d.GetListOfKeys():
      c = key.ReadObj()
      if not isinstance(c, ROOT.TCanvas):
        raise TypeError("Unknown object {} in {}/controlPlots".format(key.GetName(), templatesfile.templatesfile()))
      c.GetListOfPrimitives()[1].SetLineColor(2)  #put the line back to red, importing style breaks this for some reason
      axisnumber = int(c.GetName().split("_")[-2].replace("projAxis", ""))
      hstack = ROOT.THStack()
      hstack.Add(c.GetListOfPrimitives()[0], "P E X0")
      hstack.Add(c.GetListOfPrimitives()[1], "hist")
      c1 = ROOT.TCanvas(c.GetName(), c.GetTitle())
      hstack.Draw("nostack")
      hstack.GetXaxis().SetTitle(axistitles[axisnumber])
      for ext in exts:
        c1.SaveAs("{}/{}.{}".format(saveasdir, key.ReadObj().GetName(), ext))

if __name__ == "__main__":
    for production in productions:
        for channel in channels:
            for analysis in analyses:
                for systematic in treesystematics:
#                    if production != "160225" or channel != "2e2mu" or analysis != "fL1" or systematic != "": continue
#                    continue
                    makecontrolplots(channel, "signal", systematic, analysis, production)
#                continue
#                if channel != "4mu" or production != "160725" or analysis != "fL1": continue
                makecontrolplots(channel, "bkg", analysis, production)
