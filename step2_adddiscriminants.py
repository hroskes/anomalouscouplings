#!/usr/bin/env python

import argparse
if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--submitjobs", action="store_true")
    p.add_argument("--filter", type=eval, default=None)
    args = p.parse_args()
    if args.filter and args.submitjobs:
        p.error("Can't have a filter for submitjobs")

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
  with cdtemp_slurm(), KeepWhileOpenFile(newfilename+".tmp") as kwof, LSF_creating(newfilename, ignorefailure=True) as LSF:
    if not kwof:
        return
    if os.path.exists(newfilename):
        return

    if sample.copyfromothersample is not None:
      if sample.CJLSTfile() != sample.copyfromothersample.CJLSTfile() and xrd.exists(sample.CJLSTfile()):
        raise ValueError("{} exists, why not use it?".format(sample.CJLSTfile()))
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
            if treewrapper.alternateweightxsecstree is not None:
                treewrapper.alternateweightxsecstree.SetDirectory(newf)
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
                    try:
                        discriminants[discriminant][0] = getattr(treewrapper, discriminant)()
                    except:
                        print "Error while calculating", discriminant
                        raise
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
        if not os.path.exists(sample.withdiscriminantsfile()) and KeepWhileOpenFile(sample.withdiscriminantsfile()+".tmp").wouldbevalid and not sample.copyfromothersample:
            njobs += 1
    for i in range(njobs):
        submitjobkwargs = {"jobname": str(i), "jobtime": "1-0:0:0"}
        if config.host == "MARCC":
            submitjobkwargs["queue"] = "lrgmem"
            submitjobkwargs["memory"] = "12000M"
        submitjob("unbuffer "+os.path.join(config.repositorydir, "step2_adddiscriminants.py"), **submitjobkwargs)

if __name__ == '__main__':
    if args.submitjobs:
        submitjobs()
    else:
        try:
            for sample in allsamples():
                 if sample.productionmode == "ggZZ" and sample.flavor == "4tau" and not sample.copyfromothersample:
                     adddiscriminants(sample)
            for sample in allsamples():
                if sample.copyfromothersample: continue
                if args.filter and not args.filter(sample): continue
                adddiscriminants(sample)
        finally:
            if config.LHE and not any(KeepWhileOpenFile(sample.withdiscriminantsfile()+".tmp").wouldbevalid for sample in allsamples()):
                deletemelastuff()
