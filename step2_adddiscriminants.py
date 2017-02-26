#!/usr/bin/env python
from array import array
from collections import OrderedDict
from helperstuff import config
from helperstuff import xrd
from helperstuff.enums import hffhypotheses, ProductionMode, productions, pythiasystematics
from helperstuff.samples import allsamples, Sample
from helperstuff.submitjob import submitjob
from helperstuff.treewrapper import TreeWrapper
from helperstuff.utilities import KeepWhileOpenFile, LSB_JOBID, LSF_creating
import os
import ROOT
import sys

definitelyexists = Sample("VBF", "0+", config.productionsforcombine[0])
if not xrd.exists(definitelyexists.CJLSTfile()):
    raise ValueError("{} does not exist!".format(definitelyexists.CJLSTfile()))

def adddiscriminants(*args):

  sample = Sample(*args)
  reweightingsamples = sample.reweightingsamples()

  filename = sample.CJLSTfile()
  newfilename = sample.withdiscriminantsfile()
  print newfilename
  with KeepWhileOpenFile(newfilename+".tmp", message=LSB_JOBID()) as kwof, LSF_creating(newfilename, ignorefailure=True) as LSF:
    if not kwof:
        return
    if os.path.exists(newfilename):
        return

    if sample.copyfromothersample is not None:
      if xrd.exists(sample.CJLSTfile()):
        raise ValueError("{} exists, why not use it?".format(sample.CJLSTfile()))
      os.symlink(sample.copyfromothersample.withdiscriminantsfile(), sample.withdiscriminantsfile())
      return

    isdummy = False
    if not xrd.exists(filename):
        isdummy = True
    else:
        f = ROOT.TFile.Open(filename)
        if not f.Get("{}/candTree".format(sample.TDirectoryname())):
            isdummy = True

    if isdummy:
        print "{} does not exist or is bad, using {}".format(filename, definitelyexists.CJLSTfile())
        #give it a tree so that it can get the format, but not fill any entries
        filename = definitelyexists.CJLSTfile()
        f = ROOT.TFile.Open(filename)

    Counters = f.Get("{}/Counters".format(sample.TDirectoryname()))
    if not Counters:
        raise ValueError("No Counters in file "+filename)

    t = ROOT.TChain("{}/candTree".format(sample.TDirectoryname()))
    t.Add(filename)

    if f.Get("{}/candTree_failed".format(sample.TDirectoryname())):
        t_failed = ROOT.TChain("{}/candTree_failed".format(sample.TDirectoryname()))
        t_failed.Add(filename)
    else:
        t_failed = None

    treewrapper = TreeWrapper(t, sample, Counters=Counters, failedtree=t_failed, isdummy=isdummy)

    if os.path.exists(newfilename):
        return

    failed = False
    try:
        newf = ROOT.TFile.Open(LSF.basename(newfilename), "recreate")
        newt = t.CloneTree(0)
        if treewrapper.effectiveentriestree is not None:
            treewrapper.effectiveentriestree.SetDirectory(newf)

        discriminants = OrderedDict()
        for discriminant in treewrapper.toaddtotree:
            discriminants[discriminant] = array('d', [0])
            newt.Branch(discriminant, discriminants[discriminant], discriminant + "/D")
        for discriminant in treewrapper.toaddtotree_int:
            discriminants[discriminant] = array('i', [0])
            newt.Branch(discriminant, discriminants[discriminant], discriminant + "/I")

        try:
            for entry in treewrapper:
                for discriminant in discriminants:
                    discriminants[discriminant][0] = getattr(treewrapper, discriminant)()
                newt.Fill()
        except:
            treewrapper.Show()
            raise
    except:
        failed = True
        raise
    finally:
        try:
            newf.Write()
            newf.Close()
        except:
            failed = True
        if failed:
            try:
                os.remove(newfilename)
            except:
                pass

def submitjobs():
    njobs = 0
    for sample in allsamples():
        if not os.path.exists(sample.withdiscriminantsfile()) and not os.path.exists(sample.withdiscriminantsfile()+".tmp"):
            njobs += 1
    for i in range(njobs):
        submitjob(os.path.join(config.repositorydir, "step2_adddiscriminants.py"), jobname=str(i), jobtime="1-0:0:0")

if __name__ == '__main__':
    if sys.argv[1:]:
        if sys.argv[1].lower() == "submitjobs":
            submitjobs(*sys.argv[2:])
        else:
            raise ValueError("Can only run '{0}' with no arguments or '{0} submitjobs'".format(sys.argv[0]))
    else:
        for production in productions:
            adddiscriminants("ggZZ", "4tau", production)  #to catch bugs early
        for sample in allsamples():
            adddiscriminants(sample)
