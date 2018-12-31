#!/usr/bin/env python

import argparse
from collections import namedtuple

class stringandlambda(namedtuple("stringandlambda", "string function")):
  def __new__(cls, string):
    return super(stringandlambda, cls).__new__(cls, string, eval(string))

if __name__ == "__main__":
  p = argparse.ArgumentParser()
  p.add_argument("--submitjobs", action="store_true")
  p.add_argument("--filter", type=stringandlambda, default=None)
  args = p.parse_args()

from array import array
from collections import OrderedDict
from helperstuff import config
from helperstuff import xrd
from helperstuff.enums import hffhypotheses, ProductionMode, productions, pythiasystematics
from helperstuff.samples import allsamples, Sample
from helperstuff.submitjob import submitjob
from helperstuff.treewrapper import TreeWrapper, TreeWrapperFactory
from helperstuff.utilities import cdtemp_slurm, deletemelastuff, KeepWhileOpenFile, LSB_JOBID, LSF_creating, mkdir_p
import os
import pipes
import ROOT
import sys

def adddiscriminants(*args):

  sample = Sample(*args)
  reweightingsamples = sample.reweightingsamples()

  newfilename = sample.withdiscriminantsfile()
  print newfilename
  mkdir_p(os.path.dirname(newfilename))

  if sample.copyfromothersample is not None:
    if sample.CJLSTfile() != sample.copyfromothersample.CJLSTfile() and xrd.exists(sample.CJLSTfile()):
      raise ValueError("{} exists, why not use it?".format(sample.CJLSTfile()))
    return

  inputfiles = []
  if xrd.exists(sample.CJLSTfile()): inputfiles.append(sample.CJLSTfile())

  with cdtemp_slurm(), KeepWhileOpenFile(newfilename+".tmp") as kwof, LSF_creating(newfilename, ignorefailure=True, inputfiles=inputfiles) as LSF:
    if not kwof:
      return
    if os.path.exists(newfilename):
      return

    treewrapper = TreeWrapperFactory(sample, LSF=LSF)

    if os.path.exists(newfilename):
      return

    failed = False
    try:
      newf = ROOT.TFile.Open(LSF.basename(newfilename), "recreate")
      if isinstance(treewrapper, TreeWrapper):
        if sample.production == "180721_2016":
          treewrapper.tree.SetBranchStatus("LHEweight_PDFVariation_*", 0)
        newt = treewrapper.tree.CloneTree(0)
        if sample.production == "180721_2016":
          treewrapper.tree.SetBranchStatus("LHEweight_PDFVariation_*", 1)
          pdfup = array('d', [0])
          pdfdn = array('d', [0])
          newt.Branch("LHEweight_PDFVariation_Up", pdfup, "LHEweight_PDFVariation_Up/D")
          newt.Branch("LHEweight_PDFVariation_Dn", pdfdn, "LHEweight_PDFVariation_Dn/D")
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
        for i, entry in enumerate(treewrapper, start=1):
          for discriminant in discriminants:
            try:
              discriminants[discriminant][0] = getattr(treewrapper, discriminant)()
            except:
              print "Error while calculating", discriminant
              raise
          if sample.production == "180721_2016":
            pdfup[0] = treewrapper.tree.LHEweight_PDFVariation_Up
            pdfdn[0] = treewrapper.tree.LHEweight_PDFVariation_Dn
          newt.Fill()
          if i % 50000 == 0 and LSB_JOBID():
            with open(LSF.basename(newfilename)): pass #access it, hopefully preventing tmp from being deleted
            with open("touch.txt", "w") as f:     pass #touch another file, same idea
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

def submitjobs(filter=stringandlambda("lambda sample: True")):
  njobs = 0
  for sample in allsamples():
    if not filter.function(sample): continue
    if not os.path.exists(sample.withdiscriminantsfile()) and KeepWhileOpenFile(sample.withdiscriminantsfile()+".tmp").wouldbevalid and not sample.copyfromothersample:
      njobs += 1
  for i in range(njobs):
    submitjobkwargs = {"jobname": str(i), "jobtime": "1-0:0:0"}
    if config.host == "MARCC":
      submitjobkwargs["queue"] = "lrgmem"
      submitjobkwargs["memory"] = "12000M"
      submitjobkwargs["docd"] = True  #since cdtemp_slurm happens in thisfile
    submitjob("unbuffer "+os.path.join(config.repositorydir, "step2_adddiscriminants.py")+" --filter "+pipes.quote(filter.string), **submitjobkwargs)

if __name__ == '__main__':
  if args.submitjobs:
    submitjobs(filter=args.filter)
  else:
    try:
      for sample in allsamples():
        if sample.productionmode == "ggZZ" and sample.flavor == "4tau" and not sample.copyfromothersample:
          adddiscriminants(sample)
      for sample in allsamples():
        if sample.copyfromothersample: continue
        if args.filter and not args.filter.function(sample): continue
        adddiscriminants(sample)
    finally:
      if not any(KeepWhileOpenFile(sample.withdiscriminantsfile()+".tmp").wouldbevalid for sample in allsamples()):
        deletemelastuff()
