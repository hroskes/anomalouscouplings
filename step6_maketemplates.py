from helperstuff import config
from helperstuff.enums import DataTree, datatrees, TemplatesFile, templatesfiles
from helperstuff.samples import Sample
import os
import subprocess
import ROOT
from time import sleep

cmssw = [int(i) for i in os.environ["CMSSW_VERSION"].split("_")[1:]]
if cmssw[0] == 8:
    raise ValueError("TemplateBuilder does not seem to work in CMSSW_8_X; the templates end up filled with NaNs.  Try CMSSW_7_4_X or CMSSW_7_6_X, I have tested that it works there.")

def touch(filename):
    open(filename, 'a').close()

def buildtemplates(*args):
    templatesfile = TemplatesFile(*args)
    print templatesfile
    if os.path.exists(templatesfile.templatesfile()):
        return
    touch(templatesfile.templatesfile())
    subprocess.call([os.path.join(config.repositorydir, "TemplateBuilder/buildTemplate.exe"), templatesfile.jsonfile()])
    if not os.path.exists(templatesfile.templatesfile()):
        raise RuntimeError("Something is wrong!  {} was not created.".format(templatesfile.templatesfile()))

def copydata(*args):
    datatree = DataTree(*args)
    f = ROOT.TFile(datatree.originaltreefile)
    t = f.candTree
    newfilename = datatree.treefile
    if os.path.exists(newfilename): f.Close(); return
    newf = ROOT.TFile(newfilename, "recreate")
    newt = t.CloneTree(0)
    for entry in t:
        if abs(t.Z1Flav * t.Z2Flav) == datatree.channel.ZZFlav and config.m4lmin < t.ZZMass < config.m4lmax and config.unblindscans:
            newt.Fill()
    newf.Write()
    f.Close()
    newf.Close()

if __name__ == "__main__":
    for templatesfile in templatesfiles:
        buildtemplates(templatesfile)
        #and copy data
    for datatree in datatrees:
        copydata(datatree)
