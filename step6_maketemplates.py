#!/usr/bin/env python

import argparse
if __name__ == "__main__":
  p = argparse.ArgumentParser()
  p.add_argument("--submitjobs", action="store_true")
  p.add_argument("--jsontoo", type=int)
  p.add_argument("--removefiles", nargs="*", default=())
  p.add_argument("--waitids", nargs="*", type=int, default=())
  p.add_argument("--filter", type=eval, default=lambda template: True)
  args = p.parse_args()
  if args.jsontoo and not args.submitjobs:
    raise ValueError("--jsontoo doesn't make sense without --submitjobs")
  if args.removefiles and not args.submitjobs:
    raise ValueError("--removefiles doesn't make sense without --submitjobs")
  if args.waitids and not args.submitjobs:
    raise ValueError("--waitids doesn't make sense without --submitjobs")

from array import array
import os
import ROOT
import subprocess
import sys
from time import sleep

from helperstuff import config
from helperstuff.discriminants import discriminants
from helperstuff.enums import TemplateGroup
from helperstuff.samples import Sample
from helperstuff.submitjob import submitjob
from helperstuff.templates import DataTree, datatrees, TemplatesFile, templatesfiles
from helperstuff.utilities import cd, KeepWhileOpenFile, KeepWhileOpenFiles, LSB_JOBID, mkdir_p

def buildtemplates(*args):
  if len(args) == 1 and not isinstance(args[0], TemplatesFile):
    tg = TemplateGroup(args[0])
    tfs = [tf for tf in templatesfiles if tf.templategroup == tg and tf.usenewtemplatebuilder]
    if not tfs: return
    with KeepWhileOpenFiles(*(_.templatesfile()+".tmp" for _ in tfs)) as kwofs:
      if not all(kwofs): return
      subprocess.check_call(["buildTemplates.py", "--use-existing-templates"] + [_.jsonfile() for _ in tfs])
      return
    return

  templatesfile = TemplatesFile(*args)
  if "Ulascan" in str(templatesfile.production): return
  if templatesfile.usenewtemplatebuilder and templatesfile.templategroup in ("background", "DATA"):
    return buildtemplates(templatesfile.templategroup)
  print templatesfile
  if templatesfile.copyfromothertemplatesfile is not None: return
  with KeepWhileOpenFile(templatesfile.templatesfile() + ".tmp") as f:
    scriptname = ["buildTemplate.exe"]
    if templatesfile.usenewtemplatebuilder:
      scriptname = ["buildTemplates.py", "--use-existing-templates"]
    if f:
      if templatesfile.usenewtemplatebuilder or not os.path.exists(templatesfile.templatesfile()):
        if not os.path.exists(templatesfile.templatesfile(firststep=True)):
          mkdir_p(os.path.dirname(templatesfile.templatesfile(firststep=True)))
          try:
            subprocess.check_call(scriptname + [templatesfile.jsonfile()])
          except:
            try:
              raise
            finally:
              if not templatesfile.usenewtemplatebuilder:
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
  f = ROOT.TFile(datatree.originaltreefile)
  t = f.candTree

  discriminants_forcerange = {d: array('d', [0]) for d in discriminants.values() if hasattr(t, d.name)}
  epsilon = float("inf")
  for (dname, dtitle, dbins, dmin, dmax, didentifier, dformula), branchaddress in discriminants_forcerange.iteritems():
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

  torun = set()
  for templatesfile in templatesfiles:
    if templatesfile.copyfromothertemplatesfile is not None: continue
    if templatesfile.templatesfile() in remove:
      remove[templatesfile.templatesfile()] = True
      remove[templatesfile.templatesfile(firststep=True)] = True
      torun.add(templatesfile)
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
      torun.add(templatesfile)

  for tf in frozenset(torun):
    if tf.usenewtemplatebuilder and tf.templategroup in ("background", "DATA"):
      torun.remove(tf)
      torun.add(tf.templategroup)
  njobs = len(torun)

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
      waitids = list(step4_makejson.submitjobs(jsontoo))
    else:
      waitids = []
    waitids += list(args.waitids)
    for i in range(njobs):
      submitjob("unbuffer "+os.path.join(config.repositorydir, "step6_maketemplates.py"), jobname=str(i), jobtime="1-0:0:0", docd=True, waitids=waitids)

if __name__ == "__main__":
  if args.submitjobs:
    submitjobs(args.removefiles, args.jsontoo)
  else:
    for templatesfile in templatesfiles:
      if not args.filter(templatesfile): continue
      with cd(config.repositorydir):
        buildtemplates(templatesfile)
      #and copy data
    for datatree in datatrees:
      copydata(datatree)
