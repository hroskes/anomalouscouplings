#!/usr/bin/env python
from helperstuff import config, style  #style needed to remove the statbox and title and adjust the axis label size
from helperstuff.templates import TemplatesFile, templatesfiles
import os
import ROOT

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

    axistitles = [_.name for _ in templatesfile.discriminants]

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
      name = (key.ReadObj().GetName()
                                     .replace("afterFill", "after1Fill")
                                     .replace("afterSmooth", "after2Smooth")
                                     .replace("afterReweight", "after3Reweight")
                                     .replace("afterFloor", "after4Floor")
                                     .replace("afterNormalization", "after5Normalization")
                                     .replace("afterMirror", "after6Mirror")
             )
      for ext in exts:
        c1.SaveAs("{}/{}.{}".format(saveasdir, name, ext))

if __name__ == "__main__":
    def thetemplatesfiles():
#        yield TemplatesFile("ggh", "VBFtagged", "2e2mu", "fa3")
#        return
        for templatesfile in templatesfiles:
             yield templatesfile
    length = len(list(thetemplatesfiles()))
    for i, templatesfile in enumerate(thetemplatesfiles(), start=1):
        makecontrolplots(templatesfile)
        print i, "/", length
