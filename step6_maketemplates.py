#!/usr/bin/env python

from array import array
import os
import ROOT
import subprocess
import sys
from time import sleep

from helperstuff import config
from helperstuff.discriminants import discriminants
from helperstuff.samples import Sample
from helperstuff.submitjob import submitjob
from helperstuff.templates import DataTree, datatrees, TemplatesFile, templatesfiles
from helperstuff.utilities import cd, KeepWhileOpenFile, LSB_JOBID

cmssw = [int(i) for i in os.environ["CMSSW_VERSION"].split("_")[1:]]
if cmssw[0] == 8:
    raise ValueError("TemplateBuilder does not seem to work in CMSSW_8_X; the templates end up filled with NaNs.  Try CMSSW_7_4_X or CMSSW_7_6_X, I have tested that it works there.")

def buildtemplates(*args):
    templatesfile = TemplatesFile(*args)
    print templatesfile
    with KeepWhileOpenFile(templatesfile.templatesfile() + ".tmp", message=LSB_JOBID()) as f:
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

def submitjobs(*args):
    remove = {}
    keepjson = False
    for arg in args:
        if arg == "keepjson":
            keepjson = True
            continue
        if not arg.endswith(".root"): arg += ".root"
        arg = os.path.basename(arg)
        filename = os.path.join(config.repositorydir, "step7_templates", arg)
        if not os.path.exists(filename):
            raise ValueError("{} does not exist!".format(filename))
        remove[filename] = False

    njobs = 0
    for templatesfile in templatesfiles:
        if templatesfile.templatesfile() in remove:
            remove[templatesfile.templatesfile()] = True
            remove[templatesfile.templatesfile(firststep=True)] = True
            njobs += 1
        elif os.path.exists(templatesfile.templatesfile()) or os.path.exists(templatesfile.templatesfile() + ".tmp"):
            pass
        else:
            njobs += 1
    if not njobs: return
    for filename, found in remove.iteritems():
        if not found:
            raise IOError("{} is not a templatesfile!".format(filename))
    with cd(config.repositorydir):
        if not keepjson:
            subprocess.check_call(["./step4_makejson.py"])
        for filename in remove:
            if os.path.exists(filename):
                os.remove(filename)
        for i in range(njobs):
            submitjob(os.path.join(config.repositorydir, "step6_maketemplates.py"), jobname=str(i), jobtime="1-0:0:0")

if __name__ == "__main__":
    if sys.argv[1:]:
        if sys.argv[1].lower() == "submitjobs":
            submitjobs(*sys.argv[2:])
        else:
            raise ValueError("Can only run '{0}' with no arguments or '{0} submitjobs'".format(sys.argv[0]))
    else:
        for templatesfile in templatesfiles:
            buildtemplates(templatesfile)
            #and copy data
        for datatree in datatrees:
            copydata(datatree)
