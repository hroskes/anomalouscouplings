#!/usr/bin/env python
import glob
import os
import sys

import ROOT

from helperstuff import config, style  #style needed to remove the statbox and title and adjust the axis label size
from helperstuff.copyplots import copyplots
from helperstuff.templates import TemplatesFile, templatesfiles

def makecontrolplots(*args, **kwargs):
    templatesfile = TemplatesFile(*args)
    f = ROOT.TFile.Open(templatesfile.templatesfile(**kwargs))
    d = f.controlPlots

    saveasdir = templatesfile.controlplotsdir(**kwargs)
    try:
        os.makedirs(saveasdir)
    except OSError:
        pass

    exts = "png", "eps", "root", "pdf"

    for ext in exts:
        for filename in glob.glob(os.path.join(saveasdir, "*."+ext)):
            os.remove(filename)

    axistitles = [_.name for _ in templatesfile.discriminants]

    for key in d.GetListOfKeys():
      c = key.ReadObj()
      if not isinstance(c, ROOT.TCanvas):
        raise TypeError("Unknown object {} in {}/controlPlots".format(key.GetName(), templatesfile.templatesfile(**kwargs)))
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
                                     .replace("afterCustomSmoothing", "after7CustomSmoothing")
             )
      for ext in exts:
        c1.SaveAs("{}/{}.{}".format(saveasdir, name, ext))

if __name__ == "__main__":
    def thetemplatesfiles():
        for templatesfile in templatesfiles:
            yield templatesfile
    length = len(list(thetemplatesfiles()))

    if len(sys.argv) == 1:
        iteration = None
    elif len(sys.argv) == 2:
        iteration = sys.argv[1]
    else:
        raise TypeError("Too many command line arguments ({})".format(len(sys.argv)))

    for i, templatesfile in enumerate(thetemplatesfiles(), start=1):
        makecontrolplots(templatesfile, iteration=iteration)
        print i, "/", length

    copyplots(os.path.join("templateprojections", "controlplots"))
