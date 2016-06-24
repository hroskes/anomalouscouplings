from helperstuff import config
from helperstuff.enums import channels, treesystematics, TemplatesFile
import os
import ROOT
import subprocess

def makecontrolplots(*args):
    templatesfile = TemplatesFile(*args)
    f = ROOT.TFile.Open(templatesfile.templatesfile())
    d = f.controlPlots

    split = os.path.split(templatesfile.templatesfile())
    saveasdir = os.path.join(config.plotsbasedir, "templateprojections", "controlplots", split[1].replace("_fa3Adap_new.root", ""))
    try:
        os.makedirs(saveasdir)
    except OSError:
        pass

    exts = "png", "eps", "root", "pdf"

    for key in d.GetListOfKeys():
      for ext in exts:
        key.ReadObj().SaveAs("{}/{}.{}".format(saveasdir, key.ReadObj().GetName(), ext))

if __name__ == "__main__":
    for channel in channels:
        for systematic in treesystematics:
            makecontrolplots(channel, "signal", systematic)
        makecontrolplots(channel, "bkg")
