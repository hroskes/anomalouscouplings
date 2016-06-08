from helperstuff.enums import Channel
import os
import ROOT
import subprocess

def makecontrolplots(flavor):
    flavor = Channel(flavor)
    f = ROOT.TFile.Open(flavor.templatesfile())
    d = f.controlPlots

    split = os.path.split(flavor.templatesfile())
    saveasdir = os.path.join(split[0], "controlplots_"+split[1].replace("_fa3Adap_new.root", ""))
    try:
        os.mkdir(saveasdir)
    except OSError:
        pass

    exts = "png", "eps", "root", "pdf"

    for key in d.GetListOfKeys():
      for ext in exts:
        key.ReadObj().SaveAs("{}/{}.{}".format(saveasdir, key.ReadObj().GetName(), ext))

if __name__ == "__main__":
    makecontrolplots("4e")
    makecontrolplots("4mu")
    makecontrolplots("2e2mu")
