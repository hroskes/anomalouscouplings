#!/usr/bin/env python

import argparse
if __name__ == "__main__":
  p = argparse.ArgumentParser()
  p.add_argument("--submitjobs", action="store_true")
  p.add_argument("--jsontoo", action="store_true")
  p.add_argument("--removefiles", nargs="*", default=())
  args = p.parse_args()
  if args.jsontoo and not args.submitjobs:
    raise ValueError("--jsontoo doesn't make sense without --submitjobs")
  if args.removefiles and not args.submitjobs:
    raise ValueError("--removefiles doesn't make sense without --submitjobs")

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
from helperstuff.utilities import cd, KeepWhileOpenFile, LSB_JOBID, mkdir_p

#cmssw = [int(i) for i in os.environ["CMSSW_VERSION"].split("_")[1:]]
#if cmssw[0] == 8:
#  raise ValueError("TemplateBuilder does not seem to work in CMSSW_8_X; the templates end up filled with NaNs.  Try CMSSW_7_4_X or CMSSW_7_6_X, I have tested that it works there.")

def buildtemplates(*args):
  templatesfile = TemplatesFile(*args)
  if "Ulascan" in str(templatesfile.production): return
  print templatesfile
  if templatesfile.copyfromothertemplatesfile is not None: return
  with KeepWhileOpenFile(templatesfile.templatesfile() + ".tmp") as f:
    if f:
      if not os.path.exists(templatesfile.templatesfile()):
        if not os.path.exists(templatesfile.templatesfile(firststep=True)):
          mkdir_p(os.path.dirname(templatesfile.templatesfile(firststep=True)))
          try:
            if templatesfile.usenewtemplatebuilder: executable = "buildTemplates.py"
            else: executable = "buildTemplate.exe"
            subprocess.call([executable, templatesfile.jsonfile()])
          except:
            try:
              raise
            finally:
              try:
                os.remove(templatesfile.templatesfile(firststep=templatesfile.hascustomsmoothing))
              except:
                pass
      if (
        os.path.exists(templatesfile.templatesfile(firststep=True))
         and not os.path.exists(templatesfile.templatesfile())
         and templatesfile.hascustomsmoothing
         ):
        try:
          bad = False
          templatesfile.docustomsmoothing()
        except:
          bad = True
          raise
        finally:
          if bad and os.path.exists(templatesfile.templatesfile(firststep=True)):  #the second part has to be true.  just to make 100% sure.
            try:
              os.remove(templatesfile.templatesfile())
            except:
              pass

      if not os.path.exists(templatesfile.templatesfile()):
        raise RuntimeError("Something is wrong!  {} was not created.".format(templatesfile.templatesfile()))

def copydata(*args):
  if len(args) == 1 and isinstance(args[0], DataTree):
    datatree = args[0]
  else:
    datatree = DataTree(*args)
  if "Ulascan" in str(datatree.production): return
  print datatree
  f = ROOT.TFile(datatree.originaltreefile)
  t = f.candTree

  discriminants_forcerange = {d: array('d', [0]) for d in discriminants.values() if hasattr(t, d.name)}
  epsilon = float("inf")
  for (dname, dtitle, dbins, dmin, dmax, didentifier), branchaddress in discriminants_forcerange.iteritems():
    t.SetBranchAddress(dname, branchaddress)
    epsilon = min(epsilon, (dmax-dmin)/dbins/1000)

  newfilename = datatree.treefile
  if os.path.exists(newfilename): f.Close(); return
  mkdir_p(os.path.dirname(newfilename))

  newf = ROOT.TFile(newfilename, "recreate")
  newt = t.CloneTree(0)
  for entry in t:
    for (dname, dtitle, dbins, dmin, dmax, didentifier), branchaddress in discriminants_forcerange.iteritems():
      branchaddress[0] = min(branchaddress[0], dmax-epsilon)
      branchaddress[0] = max(branchaddress[0], dmin)
      assert dmin <= branchaddress[0] <= dmax-epsilon
    if datatree.passescut(t):
      newt.Fill()
  print newt.GetEntries()
  newf.Write()
  f.Close()
  newf.Close()

def submitjobs(removefiles, jsontoo=False):
  remove = {}
  for filename in removefiles:
    if not filename.endswith(".root"): filename += ".root"
    filename = os.path.basename(filename)
    filename = os.path.join(config.repositorydir, "step7_templates", filename)
    if not os.path.exists(filename):
      raise ValueError("{} does not exist!".format(filename))
    remove[filename] = False
    if jsontoo:
      jsonfilename = os.path.relpath(filename, config.repositorydir).replace("step7_templates", "step5_json").replace(".root", ".json")
      if os.path.exists(jsonfilename):
        remove[jsonfilename] = False

  njobs = 0
  for templatesfile in templatesfiles:
    if templatesfile.copyfromothertemplatesfile is not None: continue
    if templatesfile.templatesfile() in remove:
      remove[templatesfile.templatesfile()] = True
      remove[templatesfile.templatesfile(firststep=True)] = True
      njobs += 1
      if jsontoo:
        if templatesfile.jsonfile() in remove:
          remove[templatesfile.jsonfile()] = True
      else:
        if not os.path.exists(templatesfile.jsonfile()):
          raise ValueError(templatesfile.jsonfile()+" doesn't exist!  Try --jsontoo.")
    elif os.path.exists(templatesfile.templatesfile()) or not KeepWhileOpenFile(templatesfile.templatesfile() + ".tmp").wouldbevalid:
      pass
    else:
      if not jsontoo:
        if not os.path.exists(templatesfile.jsonfile()):
          raise ValueError(templatesfile.jsonfile()+" doesn't exist!  Try --jsontoo.")
      njobs += 1
  if not njobs: return
  for filename, found in remove.iteritems():
    if not found:
      raise IOError("{} is not a templatesfile!".format(filename))

  with cd(config.repositorydir):
    for filename in remove:
      if os.path.exists(filename):
        os.remove(filename)
    if jsontoo:
      sys.dont_write_bytecode = True
      import step4_makejson
      waitids = list(step4_makejson.submitjobs(5))
    else:
      waitids = []
    for i in range(njobs):
      submitjob("unbuffer "+os.path.join(config.repositorydir, "step6_maketemplates.py"), jobname=str(i), jobtime="1-0:0:0", docd=True, waitids=waitids)

if __name__ == "__main__":
  if args.submitjobs:
    submitjobs(args.removefiles, args.jsontoo)
  else:
    for templatesfile in templatesfiles:
      buildtemplates(templatesfile)
      #and copy data
    for datatree in datatrees:
      copydata(datatree)
