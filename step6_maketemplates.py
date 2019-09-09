#!/usr/bin/env python

import argparse
if __name__ == "__main__":
  p = argparse.ArgumentParser()
  p.add_argument("--submitjobs", action="store_true")
  p.add_argument("--jsontoo", type=int)
  p.add_argument("--removefiles", nargs="*", default=())
  p.add_argument("--waitids", nargs="*", type=int, default=())
  p.add_argument("--filter", type=eval, default=lambda template: True)
  p.add_argument("--start-with-bin", type=int, nargs=3)
  p.add_argument("--nthreads", type=int, default=8)
  p.add_argument("--on-queue", action="store_true", help=argparse.SUPPRESS)
  args = p.parse_args()
  if args.on_queue:
    args.jsontoo = None
    args.removefiles = args.waitids = ()
    args.submitjobs = False
  if args.jsontoo and not args.submitjobs:
    raise ValueError("--jsontoo doesn't make sense without --submitjobs")
  if args.removefiles and not args.submitjobs:
    raise ValueError("--removefiles doesn't make sense without --submitjobs")
  if args.waitids and not args.submitjobs:
    raise ValueError("--waitids doesn't make sense without --submitjobs")

from array import array
import os
import pipes
import ROOT
import subprocess
import sys
from time import sleep

from helperstuff import config
from helperstuff.discriminants import discriminants
from helperstuff.enums import Production, TemplateGroup
from helperstuff.samples import Sample
from helperstuff.submitjob import submitjob
from helperstuff.templates import DataTree, datatrees, TemplatesFile, templatesfiles
from helperstuff.utilities import cd, KeepWhileOpenFile, KeepWhileOpenFiles, LSB_JOBID, mkdir_p, TFile

def buildtemplates(*args, **kwargs):
  morebuildtemplatesargs = kwargs.pop("morebuildtemplatesargs", [])
  assert not kwargs, kwargs

  if len(args) == 2:
    tg = TemplateGroup(args[0])
    production = Production(args[1])
    tfs = [tf for tf in templatesfiles if tf.templategroup == tg and tf.production == production and not (os.path.exists(tf.templatesfile()) and os.path.exists(tf.templatesfile().replace(".root", ".done")))]
    if tg in ("ggh", "vbf", "zh", "wh", "vh"):
      tfs = [tf for tf in tfs if tf.shapesystematic != ""]
    if not tfs: return
    tfs = tfs[:100]
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 5000))
    with KeepWhileOpenFiles(*(_.templatesfile()+".tmp" for _ in tfs)) as kwofs:
      if not all(kwofs): return
      subprocess.check_call(["buildTemplates.py", "--use-existing-templates"] + [_.jsonfile() for _ in tfs] + morebuildtemplatesargs)
      return
    return

  templatesfile = TemplatesFile(*args)
  if "Ulascan" in str(templatesfile.production): return
  if templatesfile.templategroup in ("background", "DATA", "tth", "bbh") or templatesfile.shapesystematic != "":
    return buildtemplates(templatesfile.templategroup, templatesfile.production)
  print templatesfile
  if templatesfile.copyfromothertemplatesfile is not None: return
  with KeepWhileOpenFile(templatesfile.templatesfile() + ".tmp") as f:
    scriptname = ["buildTemplates.py", "--use-existing-templates"] + morebuildtemplatesargs
    if f:
      if (
        not os.path.exists(templatesfile.templatesfile())
        or not os.path.exists(templatesfile.templatesfile().replace(".root", ".done"))
      ):
        if not os.path.exists(templatesfile.templatesfile(firststep=True)):
          mkdir_p(os.path.dirname(templatesfile.templatesfile(firststep=True)))
          subprocess.check_call(scriptname + [templatesfile.jsonfile()])
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

  newfilename = datatree.treefile
  with KeepWhileOpenFile(newfilename + ".tmp") as kwof:
    if not kwof: return
    if os.path.exists(newfilename): return

    with TFile(datatree.originaltreefile) as f, TFile(newfilename, "recreate", deleteifbad=True) as newf:
      t = f.candTree

      discriminants_forcerange = {d: array('d', [0]) for d in discriminants.values() if hasattr(t, d.name)}
      epsilon = float("inf")
      for (dname, dtitle, dbins, dmin, dmax, didentifier, dformula), branchaddress in discriminants_forcerange.iteritems():
        t.SetBranchAddress(dname, branchaddress)
        epsilon = min(epsilon, (dmax-dmin)/dbins/1000)

      mkdir_p(os.path.dirname(newfilename))

      newt = t.CloneTree(0)
      for entry in t:
        for (dname, dtitle, dbins, dmin, dmax, didentifier, dformula), branchaddress in discriminants_forcerange.iteritems():
          branchaddress[0] = min(branchaddress[0], dmax-epsilon)
          branchaddress[0] = max(branchaddress[0], dmin)
          assert dmin <= branchaddress[0] <= dmax-epsilon
        if datatree.passescut(t):
          newt.Fill()
      print newt.GetEntries()
      newf.Write()
      f.Close()
      newf.Close()

def submitjobs(args):
  remove = {}
  for filename in args.removefiles:
    if not filename.endswith(".root"): filename += ".root"
    filename = os.path.basename(filename)
    filename = os.path.join(config.repositorydir, "step7_templates", filename)
    if not os.path.exists(filename):
      raise ValueError("{} does not exist!".format(filename))
    remove[filename] = False
    if args.jsontoo:
      jsonfilename = os.path.relpath(filename, config.repositorydir).replace("step7_templates", "step5_json").replace(".root", ".json")
      if os.path.exists(jsonfilename):
        remove[jsonfilename] = False

  torun = set()
  for templatesfile in templatesfiles:
    if templatesfile.copyfromothertemplatesfile is not None: continue
    if not args.filter(templatesfile): continue
    if templatesfile.templatesfile() in remove:
      remove[templatesfile.templatesfile()] = True
      remove[templatesfile.templatesfile(firststep=True)] = True
      torun.add(templatesfile)
      if args.jsontoo:
        if templatesfile.jsonfile() in remove:
          remove[templatesfile.jsonfile()] = True
      else:
        if not os.path.exists(templatesfile.jsonfile()):
          raise ValueError(templatesfile.jsonfile()+" doesn't exist!  Try --jsontoo.")
    elif (
      (
        os.path.exists(templatesfile.templatesfile())
        and os.path.exists(templatesfile.templatesfile().replace(".root", ".done"))
      )
      or not KeepWhileOpenFile(templatesfile.templatesfile() + ".tmp").wouldbevalid
    ):
      pass
    else:
      if not args.jsontoo:
        if not os.path.exists(templatesfile.jsonfile()):
          raise ValueError(templatesfile.jsonfile()+" doesn't exist!  Try --jsontoo.")
      torun.add(templatesfile)

  for tf in frozenset(torun):
    if tf.templategroup in ("background", "DATA", "tth", "bbh"):
      torun.remove(tf)
      torun.add((tf.templategroup, tf.production))
    elif tf.shapesystematic != "":
      torun.remove(tf)
      torun.add((tf.templategroup, tf.production, "systematic"))
  njobs = len(torun)

  if not njobs: return
  for filename, found in remove.iteritems():
    if not found:
      raise IOError("{} is not a templatesfile!".format(filename))

  with cd(config.repositorydir):
    for filename in remove:
      if os.path.exists(filename):
        os.remove(filename)
    if args.jsontoo:
      sys.dont_write_bytecode = True
      import step4_makejson
      waitids = list(step4_makejson.submitjobs(args.jsontoo))
    else:
      waitids = []
    waitids += list(args.waitids)
    for i in range(njobs):
      submitjob("unbuffer "+os.path.join(config.repositorydir, "step6_maketemplates.py")+" --on-queue " + " ".join(pipes.quote(_) for _ in sys.argv[1:]), jobname=str(i), jobtime="2-0:0:0", docd=True, waitids=waitids, memory="{}M".format(args.nthreads*6000), nthreads=args.nthreads)

if __name__ == "__main__":
  if args.submitjobs:
    submitjobs(args)
  else:
    morebuildtemplatesargs = []
    if args.start_with_bin: morebuildtemplatesargs += ["--start-with-bin"] + [str(_) for _ in args.start_with_bin]
    morebuildtemplatesargs += ["--nthreads", str(args.nthreads)]
    for templatesfile in templatesfiles:
      if not args.filter(templatesfile): continue
      with cd(config.repositorydir):
        buildtemplates(templatesfile, morebuildtemplatesargs=morebuildtemplatesargs)
      #and copy data
#    for datatree in datatrees:
#      copydata(datatree)
