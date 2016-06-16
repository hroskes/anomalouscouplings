from helperstuff import config
from helperstuff.enums import Channel
import os
import ROOT
import subprocess

def makecontrolplots(flavor, isbkg):
    flavor = Channel(flavor)
    f = ROOT.TFile.Open(flavor.templatesfile(isbkg))
    d = f.controlPlots

    split = os.path.split(flavor.templatesfile(isbkg))
    saveasdir = os.path.join(config.plotsbasedir, "templateprojections", "controlplots_"+split[1].replace("_fa3Adap_new.root", ""))
    try:
        os.makedirs(saveasdir)
    except OSError:
        pass

    exts = "png", "eps", "root", "pdf"

    for key in d.GetListOfKeys():
      for ext in exts:
        key.ReadObj().SaveAs("{}/{}.{}".format(saveasdir, key.ReadObj().GetName(), ext))

if __name__ == "__main__":
    makecontrolplots("4e", True)
    makecontrolplots("4mu", True)
    makecontrolplots("2e2mu", True)
    makecontrolplots("4e", False)
    makecontrolplots("4mu", False)
    makecontrolplots("2e2mu", False)
