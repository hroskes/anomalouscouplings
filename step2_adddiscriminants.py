#!/usr/bin/env python
from array import array
from collections import OrderedDict
from helperstuff import config
from helperstuff import xrd
from helperstuff.enums import hffhypotheses, ProductionMode, productions, pythiasystematics
from helperstuff.samples import allsamples, Sample
from helperstuff.submitjob import submitjob
if config.LHE:
    from helperstuff.lhewrapper import LHEWrapper as TreeWrapper
else:
    from helperstuff.treewrapper import TreeWrapper
from helperstuff.utilities import cdtemp_slurm, deletemelastuff, KeepWhileOpenFile, LSB_JOBID, LSF_creating
import os
import ROOT
import sys

def adddiscriminants(*args):

  sample = Sample(*args)
  reweightingsamples = sample.reweightingsamples()

  newfilename = sample.withdiscriminantsfile()
  print newfilename
  with cdtemp_slurm(), KeepWhileOpenFile(newfilename+".tmp", message=LSB_JOBID()) as kwof, LSF_creating(newfilename, ignorefailure=True) as LSF:
    if not kwof:
        return
    if os.path.exists(newfilename):
        return

    if sample.copyfromothersample is not None:
      if xrd.exists(sample.CJLSTfile()):
        raise ValueError("{} exists, why not use it?".format(sample.CJLSTfile()))
      os.symlink(sample.copyfromothersample.withdiscriminantsfile(), sample.withdiscriminantsfile())
      return

    treewrapper = TreeWrapper(sample)

    if os.path.exists(newfilename):
        return

    failed = False
    try:
        newf = ROOT.TFile.Open(LSF.basename(newfilename), "recreate")
        if not config.LHE:
            newt = treewrapper.tree.CloneTree(0)
            if treewrapper.effectiveentriestree is not None:
                treewrapper.effectiveentriestree.SetDirectory(newf)
        else:
            newt = ROOT.TTree("candTree", "candTree")

        discriminants = OrderedDict()
        for discriminant in treewrapper.toaddtotree:
            discriminants[discriminant] = array('d', [0])
            newt.Branch(discriminant, discriminants[discriminant], discriminant + "/D")
        for discriminant in treewrapper.toaddtotree_int:
            discriminants[discriminant] = array('i', [0])
            newt.Branch(discriminant, discriminants[discriminant], discriminant + "/I")
        for discriminant in treewrapper.toaddtotree_float:
            discriminants[discriminant] = array('f', [0])
            newt.Branch(discriminant, discriminants[discriminant], discriminant + "/F")

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
        try:
            for production in productions:
                if not config.LHE:
                    adddiscriminants("ggZZ", "4tau", production)  #to catch bugs early
            for sample in allsamples():
                adddiscriminants(sample)
        finally:
            if config.LHE and not any(os.path.exists(sample.withdiscriminantsfile()+".tmp") for sample in allsamples()):
                deletemelastuff()
