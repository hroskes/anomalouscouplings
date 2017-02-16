#!/usr/bin/env python

from array import array
import os
import ROOT
import subprocess
from time import sleep

from helperstuff import config
from helperstuff.discriminants import discriminants
from helperstuff.samples import Sample
from helperstuff.templates import DataTree, datatrees, TemplatesFile, templatesfiles
from helperstuff.utilities import KeepWhileOpenFile

cmssw = [int(i) for i in os.environ["CMSSW_VERSION"].split("_")[1:]]
if cmssw[0] == 8:
    raise ValueError("TemplateBuilder does not seem to work in CMSSW_8_X; the templates end up filled with NaNs.  Try CMSSW_7_4_X or CMSSW_7_6_X, I have tested that it works there.")

def buildtemplates(*args):
    templatesfile = TemplatesFile(*args)
    print templatesfile
    with KeepWhileOpenFile(templatesfile.templatesfile() + ".tmp") as f:
        if f:
            if not os.path.exists(templatesfile.templatesfile()):
                if not os.path.exists(templatesfile.templatesfile(firststep=True)):
                    try:
                        subprocess.call([os.path.join(config.repositorydir, "TemplateBuilder/buildTemplate.exe"), templatesfile.jsonfile()])
                    except:
                        try:
                            os.remove(templatesfile.templatesfile(firststep=True))
                        except:
                            pass
                        raise
            if (
                templatesfile.hascustomsmoothing
                 and os.path.exists(templatesfile.templatesfile(firststep=True))
                 and not os.path.exists(templatesfile.templatesfile())
               ):
                try:
                    templatesfile.docustomsmoothing()
                except:
                    if os.path.exists(templatesfile.templatesfile(firststep=True)):  #this has to be true.  just to make 100% sure.
                        try:
                            os.remove(templatesfile.templatesfile())
                        except:
                            pass
                    raise

            if not os.path.exists(templatesfile.templatesfile()):
                raise RuntimeError("Something is wrong!  {} was not created.".format(templatesfile.templatesfile()))

def copydata(*args):
    if len(args) == 1 and isinstance(args[0], DataTree):
        datatree = args[0]
    else:
        datatree = DataTree(*args)
    print datatree
    f = ROOT.TFile(datatree.originaltreefile)
    t = f.candTree

    discriminants_forcerange = {d: array('d', [0]) for d in discriminants.values() if hasattr(t, d.name)}
    for (dname, dtitle, dbins, dmin, dmax), branchaddress in discriminants_forcerange.iteritems():
        t.SetBranchAddress(dname, branchaddress)
    epsilon = (dmax-dmin)/dbins/1000

    newfilename = datatree.treefile
    if os.path.exists(newfilename): f.Close(); return

    newf = ROOT.TFile(newfilename, "recreate")
    newt = t.CloneTree(0)
    for entry in t:
        for (dname, dtitle, dbins, dmin, dmax), branchaddress in discriminants_forcerange.iteritems():
            branchaddress[0] = min(branchaddress[0], dmax-epsilon)
            branchaddress[0] = max(branchaddress[0], dmin)
            assert dmin <= branchaddress[0] <= dmax-epsilon
        if datatree.passescut(t):
            newt.Fill()
    print newt.GetEntries()
    newf.Write()
    f.Close()
    newf.Close()

if __name__ == "__main__":
    buildtemplates("tth", "4mu", "170203", "Untagged", "fL1")
if False:
    for templatesfile in templatesfiles:
        buildtemplates(templatesfile)
        #and copy data
    for datatree in datatrees:
        copydata(datatree)
