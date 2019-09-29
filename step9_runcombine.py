#!/usr/bin/env python
import getpass
import glob
from itertools import izip, izip_longest, permutations, product
import json
import os
import pipes
import re
import shutil
import subprocess
import sys

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True #https://root-forum.cern.ch/t/pyroot-crashes-when-in-arguments/25379/3

from Alignment.OfflineValidation.TkAlAllInOneTool.helperFunctions import replaceByMap  #easiest place to get it

from helperstuff import config, utilities
from helperstuff.combinehelpers import Luminosity
from helperstuff.datacard import makeDCsandWSs
from helperstuff.enums import Analysis, categories, Category, Channel, channels, Production, ProductionMode
from helperstuff.plotlimits import plotlimits, plotlimits2D, plottitle
from helperstuff.submitjob import submitjob
from helperstuff.utilities import cd, deprecate, LSB_JOBID, mkdir_p, PlotCopier, requirecmsenv, tfiles

requirecmsenv(os.path.join(config.repositorydir, "CMSSW_10_2_5"))

combinecardstemplate = r"""
combineCards.py .oO[cardstocombine]Oo. > .oO[combinecardsfile]Oo.
"""

createworkspacetemplate = r"""
set -euo pipefail &&
unbuffer text2workspace.py -m 125 .oO[combinecardsfile]Oo. -P .oO[physicsmodel]Oo. \
                           .oO[physicsoptions]Oo. -o .oO[workspacefile]Oo. -v 7 .oO[turnoff]Oo. \
                           |& tee log.text2workspace.oO[workspacefileappend]Oo.
"""
runcombinetemplate = r"""
set -euo pipefail &&
.oO[combine]Oo. -M .oO[method]Oo. -d .oO[workspacefile]Oo. --robustFit=.oO[robustfit]Oo. \
                .oO[morecombineoptions]Oo. \
                --setParameterRanges .oO[parameterranges]Oo. -m 125 .oO[setPOI]Oo. \
                -n _.oO[append]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo. .oO[selectpoints]Oo. \
                .oO[saveorloadworkspace]Oo. \
                --X-rtd OPTIMIZE_BOUNDS=0 --X-rtd TMCSO_AdaptivePseudoAsimov=0 --X-rtd MINIMIZER_analytic \
                --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_MaxCalls=999999999 \
                .oO[-t -1]Oo. --setParameters .oO[setparameters]Oo. -V -v 3 --saveNLL \
                --setParametersForGrid=.oO[setparametersforgrid]Oo. \
|& tee .oO[logfile]Oo.
"""

plotimpactscommand = r"""
plotImpacts.py -i .oO[impactsfilename3]Oo. -o .oO[saveasdir]Oo./impacts_.oO[append]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo. &&
gm convert .oO[saveasdir]Oo./impacts_.oO[append]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo..pdf .oO[saveasdir]Oo./impacts_.oO[append]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo..png
"""

morecombineoptions = {
    "MultiDimFit": "--algo .oO[algo]Oo. --points .oO[internalnpoints]Oo. --floatOtherPOIs=.oO[floatotherpois]Oo. --alignEdges=1  .oO[savemu]Oo. .oO[savespecifiednuis]Oo. --saveInactivePOI=1 .oO[freezeparameters]Oo.",
    "FitDiagnostics": "",
    "Impacts": ".oO[impactsstep.oO[impactsstep]Oo.]Oo."
}

diffnuisancescommand = r"""
                  $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py \
                  -a fitDiagnostics_.oO[append]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo..root \
                  --poi .oO[POI]Oo. |& tee fitDiagnostics.oO[append]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo..txt
"""

def check_call_test(*args, **kwargs):
    print args[0]
#uncomment this for testing purposes
#subprocess.check_call = subprocess.check_output = check_call_test

def runscan(repmap, submitjobs, directory=None):
  originalrepmap = repmap

  if submitjobs and repmap["method"] == "MultiDimFit":
    mkdir_p("jobs")
    npoints = int(repmap["npoints"])
    jobids = set()
    individualfilenames = set()

    repmap_final = repmap.copy()
    repmap_final["selectpoints"] = ""
    finalfilename = replaceByMap(".oO[filename]Oo.", repmap_final)
    if os.path.exists(finalfilename): return

    repmap_initial = repmap.copy()
    repmap_initial.update({
      "npoints": "101",
      "scanrange": {
        "CMS_zz4l_fai1": "-1,1",
        "CMS_zz4l_fai2": "-1,1",
        "CMS_zz4l_fai3": "-1,1",
        "CMS_zz4l_fai4": "-1,1",
      }[repmap["POI"]],
      "scanrangeappend": ".oO[pointindex]Oo.",
      "expectfai": "0.0",
      "pointindex": "_firststep",
      "selectpoints": "--firstPoint 1 --lastPoint 0",
      "saveorloadworkspace": "--saveWorkspace",
    })

    initialjobids = []
    print replaceByMap("jobs/.oO[filename]Oo.", repmap_initial)
    if not utilities.existsandvalid(replaceByMap("jobs/.oO[filename]Oo.", repmap_initial), "w"):
      job = "{} runscan {} directory={}".format(
                                                os.path.join(config.repositorydir, "./step9_runcombine.py"),
                                                pipes.quote(json.dumps(repmap_initial)),
                                                pipes.quote(os.getcwd()),
                                               )
      submitjobkwargs = {
        "jobname": replaceByMap(".oO[expectfaiappend]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo...oO[exporobs]Oo.", repmap_initial).lstrip("_"),
        "jobtime": "1-0:0:0",
        "outputfile": replaceByMap("jobs/joblog.oO[expectfaiappend]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo...oO[exporobs]Oo.", repmap_initial),
        "morerepmap": repmap_initial
      }

      jobid = submitjob(job, **submitjobkwargs)
      initialjobids.append(jobid)

    for i in range(npoints):
      repmap_i = repmap.copy()
      repmap_i.update({
        "selectpoints": "--firstPoint .oO[pointindexplusoffset]Oo. --lastPoint .oO[pointindexplusoffset]Oo.",
        "pointindex": str(i),
        "pointindexplusoffset": str(i+int(replaceByMap(".oO[pointoffset]Oo.", repmap_i))),
        "saveorloadworkspace": "--snapshotName MultiDimFit --skipInitialFit",
        "workspacefile": replaceByMap(".oO[filename]Oo.", repmap_initial),
      })
      replaceByMap(runcombinetemplate, repmap_i) #sanity check that replaceByMap works
      individualfilenames.add(replaceByMap("jobs/.oO[filename]Oo.", repmap_i))

      ############################
      #compatibility
      if os.path.exists(replaceByMap(".oO[filename]Oo.", repmap_i)):
        shutil.move(replaceByMap(".oO[filename]Oo.", repmap_i), "jobs/")
      ############################

      if not utilities.KeepWhileOpenFile(replaceByMap("jobs/.oO[filename]Oo..tmp", repmap_i)).wouldbevalid:
        continue

      utilities.existsandvalid(replaceByMap("jobs/.oO[filename]Oo.", repmap_i), "limit")

      if os.path.exists(replaceByMap("jobs/.oO[filename]Oo.", repmap_i)):
        continue
      job = "{} runscan {} directory={}".format(
                                                os.path.join(config.repositorydir, "./step9_runcombine.py"),
                                                pipes.quote(json.dumps(repmap_i)),
                                                pipes.quote(os.getcwd()),
                                               )
      submitjobkwargs = {
        "jobname": replaceByMap(".oO[expectfaiappend]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo...oO[exporobs]Oo.", repmap_i).lstrip("_"),
        "jobtime": "1-0:0:0",
        "outputfile": replaceByMap("jobs/joblog.oO[expectfaiappend]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo...oO[exporobs]Oo.", repmap_i),
        "morerepmap": repmap_i,
        "waitids": initialjobids,
      }
      if config.host == "MARCC":
        submitjobkwargs["queue"] = "lrgmem"
        submitjobkwargs["memory"] = "10G"
      jobid = submitjob(job, **submitjobkwargs)
      jobids.add(jobid)


    haddjob = " ".join([os.path.join(config.repositorydir, "helperstuff", "hadd.py"), finalfilename] + list(individualfilenames))
    submithaddjobkwargs = {
      "jobname": "hadd",
      "jobtime": "1:0:0",
      "waitids": jobids,
      "docd": True,
    }
    if config.host == "MARCC":
      submithaddjobkwargs["queue"] = "lrgmem"
      submithaddjobkwargs["memory"] = "10G"

    haddjobid = submitjob(haddjob, **submithaddjobkwargs)
    return haddjobid

  else:
    cwd = os.path.realpath(os.getcwd())
    if directory is None: directory = cwd
    directory = os.path.realpath(directory)
    if cwd != directory:
#      if utilities.LSB_JOBID() is None:
#        raise ValueError("Should call runscan from directory except in a batch job")
      shutil.copy(
        os.path.join(
          directory,
          "jobs" if "higgsCombine" in replaceByMap(".oO[workspacefile]Oo.", repmap) else "",
          replaceByMap(".oO[workspacefile]Oo.", repmap)
        ),
        "."
      )
      if utilities.LSB_JOBID() is None:
         cdto = os.path.basename(directory.rstrip("/"))
      else:
         directory = os.path.join(directory, "jobs")
         cdto = "."
    else:
      cdto = "."

    repmap = repmap.copy()
    if "selectpoints" not in repmap:
      repmap["selectpoints"] = ""

    if repmap["method"] == "Impacts":
      tmp = repmap.copy()
      tmp["impactsstep"] = "1"

      tmp["seedforimpacts"] = ".123456"
      linkfrom = None
      if utilities.existsandvalid(replaceByMap(".oO[filename]Oo.", tmp)):
        linkfrom = replaceByMap(".oO[filename]Oo.", tmp)
      tmp["seedforimpacts"] = ""
      if linkfrom and not utilities.existsandvalid(replaceByMap(".oO[filename]Oo.", tmp)):
        os.symlink(linkfrom, replaceByMap(".oO[filename]Oo.", tmp))

      if os.path.exists(replaceByMap(".oO[filename]Oo.", tmp)):
        tmp["impactsstep"] = "2"
        nuisances = {"RV": None, "RF": None}
        with open(replaceByMap(".oO[combinecardsfile]Oo.", repmap)) as f:
          for line in f:
            try:
              if not line.strip(): continue
              if line[0] in "-#": continue
              split = line.split()
              if split[0] in "Combination imax jmax kmax shapes bin observation process rate": continue
              if split[1] in "lnN param shape1": nuisances[split[0]] = None; continue
            except:
              print "\n\n\nError when reading line\n\n"+line+"\n\n"
              raise
            raise ValueError("Don't know what to do with line:\n\n"+line)

        exist = {}
        for nuisance in nuisances.keys():
          tmp["nuisanceforimpacts"] = nuisance

          tmp["seedforimpacts"] = ".123456"
          linkfrom = None
          if utilities.existsandvalid(replaceByMap(".oO[filename]Oo.", tmp)):
            linkfrom = replaceByMap(".oO[filename]Oo.", tmp)
          tmp["seedforimpacts"] = ""
          if linkfrom and not utilities.existsandvalid(replaceByMap(".oO[filename]Oo.", tmp)):
            os.symlink(linkfrom, replaceByMap(".oO[filename]Oo.", tmp))

          nuisances[nuisance] = utilities.existsandvalid(replaceByMap(".oO[filename]Oo.", tmp), "limit")

        if all(nuisances.itervalues()):
          tmp["impactsstep"] = "3"
        elif any(nuisances.itervalues()):
          raise ValueError("Files exist for some but not all nuisances.\n\nExist:\n{}\n\nDon't exist:\n{}".format(
            "\n".join("  "+os.path.basename(k) for k, v in nuisances.iteritems() if v),
            "\n".join("  "+os.path.basename(k) for k, v in nuisances.iteritems() if not v),
          ))
        else:  #none exist, keep step 2
          repmap["nuisanceforimpacts"] = sorted(nuisances.iterkeys())[0]

      if "impactsstep" in repmap and repmap["impactsstep"] != tmp["impactsstep"]:
        raise RuntimeError("Trying to run impacts step {}, but {} does not exist".format(repmap["impactsstep"], replaceByMap(".oO[filename]Oo.", tmp)))
      repmap["impactsstep"] = tmp["impactsstep"]
      del tmp

    filename = replaceByMap(".oO[filename]Oo.", repmap)

    tmpfile = os.path.join(directory, filename+".tmp")
    logfile = replaceByMap(".oO[logfile]Oo.", repmap)

    if not os.path.exists(filename):
      with cd(cdto), \
           utilities.OneAtATime(tmpfile, 30), \
           utilities.LSF_creating(os.path.join(directory, filename), skipifexists=True), \
           utilities.LSF_creating(os.path.join(directory, logfile), skipifexists=True):
        if not os.path.exists(filename):
          subprocess.check_call(replaceByMap(runcombinetemplate, repmap), shell=True)

    if repmap["method"] == "Impacts" and repmap["impactsstep"] != "3":
      newrepmap = originalrepmap.copy()
      newrepmap["impactsstep"] = str(int(repmap["impactsstep"]) + 1)
      runscan(repmap=originalrepmap, submitjobs=submitjobs, directory=directory)


def runcombine(analysis, foldername, **kwargs):
    global ntry
    ntry += 1
    inputargs = analysis, foldername
    inputkwargs = kwargs.copy()

    usechannels = channels
    usecategories = categories
    expectvalues = [0.0]
    plotname = "limit"
    legendposition = (.2, .7, .6, .9)
    CLtextposition = "left"
    productions = config.productionsforcombine
    usesystematics = True
    runobs = config.unblindscans
    subdirectory = ""
    subdirectory2 = ""
    defaultscanrange = (101, -1.0, 1.0)
    scanranges = [defaultscanrange]
    defaultusesignalproductionmodes = usesignalproductionmodes = {ProductionMode(p) for p in ("ggH", "VBF", "ZH", "WH", "ttH", "bbH")}
    usebkg = True
    fixmuV = fixmuf = fixfai = False
    plotnuisances = []
    defaultalgo = algo = "grid"
    robustfit = False
    POI = None
    submitjobs = False
    alsocombinename = None
    alsocombine = []
    sqrts = None
    lumitype = "fordata"
    CMStext = "Preliminary"
    drawCMS = True
    scanfai = analysis
    faifor = "decay"
    xaxislimits = None
    method = "MultiDimFit"
    plotlimitskwargs = {"combinelogs": []}
    ntoys = -1
    plotcopier = None
    freeze = {}
    faiorder = None
    floatothers = True
    setparametersforgrid = None
    onlyworkspace = False
    for kw, kwarg in kwargs.iteritems():
        if kw == "channels":
            usechannels = [Channel(c) for c in kwarg.split(",")]
        elif kw == "categories":
            usecategories = [Category(c) for c in kwarg.split(",")]
        elif kw == "expectvalues":
            try:
                expectvalues = [fai if fai=="minimum" else float(fai) for fai in kwarg.split(",")]
            except ValueError:
                if kwarg == "":
                    expectvalues = []
                else:
                    raise ValueError("expectvalues has to contain floats separated by commas!")
        elif kw == "plotname":
            plotname = kwarg
        elif kw == "CLtextposition":
            CLtextposition = kwarg
        elif kw == "legendposition":
            try:
                legendposition = [float(a) for a in kwarg.split(",")]
                if len(legendposition) != 4: raise ValueError
            except ValueError:
                raise ValueError("legendposition has to contain 4 floats separated by commas!")
        elif kw == "productions":
            productions = [Production(p) for p in kwarg.split(",")]
        elif kw == "usesystematics":
            usesystematics = bool(int(kwarg))
        elif kw == "scanranges":
            scanranges = []
            kwarg = kwarg.replace(":", ";")
            for _ in kwarg.split(";"):
                if not _: continue
                try:
                    split = _.split(",")
                    scanrange = tuple([int(split[0])] + [float(a) for a in split[1:]])
                    if len(scanrange) != 3: raise ValueError
                    scanranges.append(scanrange)
                except ValueError:
                    raise ValueError("each scanrange has to contain an int and 2 floats separated by commas!")
        elif kw == "usesignalproductionmodes":
            usesignalproductionmodes = [ProductionMode(p) for p in sorted(kwarg.split(","))]
        elif kw == "subdirectory":
            subdirectory = kwarg
        elif kw == "subdirectory2":
            subdirectory2 = kwarg
        elif kw == "usebkg":
            if kwarg.lower() == "true":
                usebkg = True
            elif kwarg.lower() == "false":
                usebkg = False
            else:
                usebkg = bool(int(kwarg))
        elif kw == "luminosity":
            lumitype = float(kwarg)
        elif kw == "fixmuV":
            fixmuV = bool(int(kwarg))
        elif kw == "fixmuf":
            fixmuf = bool(int(kwarg))
        elif kw == "fixfai":
            fixfai = bool(int(kwarg))
        elif kw == "plotnuisances":
            plotnuisances += kwarg.split(",")
        elif kw == "plotmus":
            if bool(int(kwarg)):
                plotnuisances += ["muV", "muf"]
        elif kw == "algo":
            algo = kwarg
        elif kw == "robustfit":
            robustfit = bool(int(kwarg))
        elif kw == "POI":
            POI = kwarg
        elif kw == "runobs":
            runobs = bool(int(kwarg))
        elif kw == "submitjobs":
            submitjobs=bool(int(kwarg))
        elif kw == "alsocombine":
            kwarg = kwarg.split(",")
            alsocombinename = kwarg[0]
            if ".root" in alsocombinename or ".txt" in alsocombinename:
                raise ValueError("For alsocombine, the first item should be a name")

            if len(kwarg) == 1:
                raise ValueError("You didn't give anything to alsocombine")

            alsocombine = []

            tocopy = []
            for _ in kwarg[1:]:
                if _.startswith("/work-zfs") and config.host == "lxplus":
                    tocopy.append(
                        config.getconfiguration("login-node02", config.marccusername)["connect"]
                        + ":" + _.replace(".txt", "*").replace(".input.root", "*").replace(".root", "*")
                    )
            if tocopy:
                assert not LSB_JOBID()
                tmpdir = "/tmp/"+getpass.getuser()+"/{}_{}/".format(analysis, foldername)
                mkdir_p(tmpdir)
                with cd(tmpdir):
                    subprocess.check_call(["rsync", "-azvPR"] + tocopy + ["."])
                for i, _ in enumerate(kwarg):
                    if _.startswith("/work-zfs"):
                        kwarg[i] = tmpdir + _


            for _ in kwarg[1:]:
                theglob = glob.glob(_)
                if not theglob:
                    raise ValueError("{} does not exist!".format(_))
                theglob = [_.replace(".input.root", ".txt") for _ in theglob]
                theglob = list(set(theglob))
                alsocombine += theglob
            del theglob

            for _ in alsocombine:
                if ".txt" not in _:
                    raise ValueError("{} in alsocombine should end with .txt or .input.root".format(_))
                for _2 in _, _.replace(".txt", ".input.root"):
                    if not os.path.exists(_):
                        raise ValueError("{} does not exist!".format(_2))
        elif kw == "sqrts":
            try:
                sqrts = [int(_) for _ in kwarg.split(",")]
            except ValueError:
                raise ValueError("sqrts has to contain ints separated by commas!")
        elif kw == "CMStext":
            CMStext = kwarg
        elif kw == "drawCMS":
            drawCMS = bool(int(kwarg))
        elif kw == "scanfai":
            scanfai = Analysis(kwarg)
        elif kw == "faifor":
            faifor = kwarg
        elif kw == "xaxislimits":
            try:
                xaxislimits = [float(_) for _ in kwarg.split(",")]
                if len(xaxislimits) != 2:
                    raise ValueError
            except ValueError:
                raise ValueError("xaxislimits has to be 2 floats separated by commas")
        elif kw in ("killpoints", "useNLLandNLL0"):
            plotlimitskwargs[kw] = kwarg
        elif kw == "method":
            method = kwarg
        elif kw == "ntoys":
            ntoys = int(kwarg)
        elif kw == "freeze":
            for thing in kwarg.split(","):
                parameter, value = thing.split(":")
                for fai in analysis.fais:
                    if parameter == str(fai):
                        parameter = "CMS_zz4l_fai{}".format(analysis.fais.index(fai)+1)
                freeze[parameter] = value
        elif kw == "plotcopier":
            plotcopier = kwarg
        elif kw == "faiorder":
            faiorder = tuple(Analysis(_) if _ != "fa1" else _ for _ in kwarg.split(","))
        elif kw == "floatothers":
            floatothers = bool(int(kwarg))
        elif kw == "setparametersforgrid":
            setparametersforgrid = kwarg.replace(":", "=")
        elif kw == "onlyworkspace":
            onlyworkspace = bool(int(kwarg))
        else:
            raise TypeError("Unknown kwarg: {}".format(kw))

    defaultPOI = "CMS_zz4l_fai{}".format(analysis.fais.index(scanfai)+1)
    if POI is None: POI = defaultPOI

    if runobs and not config.unblindscans:
        raise TypeError("Can't unblind scans!")
    if runobs and lumitype != "fordata":
        raise TypeError("For unblindscans, if you want to adjust the luminosity do it in the Production class (in enums.py)")
    if scanfai not in analysis.fais:
        raise ValueError("scanfai for "+str(analysis)+" has to be one of " + " or ".join(str(_) for _ in list(analysis.fais)) + ", "+str(scanfai)+" is invalid")

    is2dscan = len(scanfai.fais) > 1
    if is2dscan:
      scanranges = list((points**2, min, max) for points, min, max in scanranges)

    internalscanranges = []
    for scanrange in scanranges:
      internalscanrange = list(scanrange) + [0] #0 is the offset
      internalscanranges.append(internalscanrange)
      assert not is2dscan
      if POI in ("CMS_zz4l_fai1", "CMS_zz4l_fai2", "CMS_zz4l_fai3", "CMS_zz4l_fai4"):
        hastoinclude = 0
      else:
        assert False, POI
      assert scanrange[2] > scanrange[1], (scanrange[2], scanrange[1])
      while not internalscanrange[1] < hastoinclude:
        internalscanrange[1] -= (internalscanrange[2] - internalscanrange[1])
        internalscanrange[3] += scanrange[0]-1
        internalscanrange[0] += scanrange[0]-1
        assert submitjobs
      while not internalscanrange[2] > hastoinclude:
        internalscanrange[2] += (internalscanrange[2] - internalscanrange[1])
        internalscanrange[0] += scanrange[0]-1
        assert submitjobs

    if submitjobs:
        if not utilities.inscreen():
            raise RuntimeError("submitjobs should be run from a screen session!")
        if "minimum" in expectvalues:
            raise ValueError("Can't run scan for minimum using submitjobs")
        if method != "MultiDimFit":
            raise ValueError("Can't do submitjobs for the "+method+" method")

    if "minimum" in expectvalues and method != "MultiDimFit":
        raise ValueError("Can't run scan for minimum using "+method)

    if len(productions) > 1 and lumitype != "fordata":
        raise TypeError("If there's >1 production, have to use lumitype == 'fordata'")

    years = [p.year for p in productions]
    if len(set(years)) != len(years):
        raise ValueError("Some of your productions are from the same year!")

    LHE = {(p.LHE, p.GEN) for p in productions}
    if len(LHE) != 1:
        raise ValueError("Some of your productions are for LHE and some are not!")
    LHE, GEN = LHE.pop()

    defaultfaiorder = tuple(sorted(analysis.fais, key=lambda x: x!=scanfai) + ["fa1"])
    if faiorder is None: faiorder = defaultfaiorder
    if set(faiorder) != set(defaultfaiorder):
      raise ValueError("faiorder doesn't include the right fais.\n{}\n{}".format(set(faiorder), set(defaultfaiorder)))
    if faiorder[0] != scanfai:
      raise ValueError("{} should be first in faiorder".format(scanfai))

    if set(usechannels) != {Channel("2e2mu")} and LHE:
        raise ValueError("For LHE analysis have to specify channels=2e2mu")

    combinecardsappend = "_lumi.oO[totallumi]Oo."
    workspacefileappend = ".oO[combinecardsappend]Oo."
    moreappend = ".oO[workspacefileappend]Oo."
    turnoff = []
    if fixmuV:
        assert False
        workspacefileappend += "_fixmuV"
    if fixmuf:
        assert False
        workspacefileappend += "_fixmuf"
    if set(usesignalproductionmodes) != set(defaultusesignalproductionmodes):
        workspacefileappend += "_"+",".join(str(p) for p in usesignalproductionmodes)
        disableproductionmodes = set(defaultusesignalproductionmodes) - set(usesignalproductionmodes)
        turnoff.append("--PO turnoff={}".format(",".join(p.combinename for p in disableproductionmodes)))
        if not any(_ in ("VBF", "ZH", "WH") for _ in usesignalproductionmodes):
            turnoff.append("--PO noRV")
        if not any(_ in ("ggH", "ttH", "bbH") for _ in usesignalproductionmodes):
            turnoff.append("--PO noRF")
    if set(usechannels) != set(channels):
        combinecardsappend += "_" + ",".join(sorted(str(c) for c in usechannels))
    if set(usecategories) != set(categories):
        combinecardsappend += "_" + ",".join(sorted(str(c) for c in usecategories))
    if not usebkg:
        workspacefileappend += "_nobkg"
        turnoff.append("--PO nobkg")
    if not usesystematics:
        freeze["allConstrainedNuisances"] = None
        moreappend += "_nosystematics"
    if algo != defaultalgo:
        moreappend += "_algo"+algo
    if robustfit:
        moreappend += "_robustfit"
    if POI != defaultPOI:
        if analysis.dimensions == 2: assert False
        moreappend += "_scan"+plottitle(POI)
    if fixfai:
        moreappend += "_fixfai"
    if not floatothers:
        moreappend += "_fixothers"
    if alsocombine:
        combinecardsappend += "_" + alsocombinename
        if sqrts is None:
            raise ValueError("Have to provide sqrts if you provide alsocombine!")
    if scanfai != analysis:
        workspacefileappend += "_scan{}".format(scanfai)
    for k, v in freeze.iteritems():
        if v is not None:
            moreappend += "_{}={}".format(k, v)
    if [str(_) for _ in faiorder] != [str(_) for _ in defaultfaiorder]:
        workspacefileappend += "_"+",".join(str(_) for _ in faiorder)
    if setparametersforgrid:
        moreappend += "_"+setparametersforgrid

    if not analysis.useboosted:
        usecategories = [c for c in usecategories if c != "Boosted"]
    if not analysis.usemorecategories:
        usecategories = [c for c in usecategories if c not in ("VBF1jtagged", "VHLepttagged")]
    if analysis.isdecayonly:
        usecategories = [c for c in usecategories if c == "Untagged"]

    if len(analysis.fais) > 1: assert not fixfai

    if sqrts is None:
        sqrts = [13]

    analysis = Analysis(analysis)
    fullfoldername = "{}_{}".format(analysis, foldername)
    totallumi = sum(float(Luminosity(p, lumitype)) for p in productions)
    saveasdir = os.path.join(config.plotsbasedir, "limits", subdirectory, fullfoldername, subdirectory2)
    try:
        os.makedirs(saveasdir)
    except OSError:
        pass

    parameterranges = {
      defaultPOI: "-1,1",
      POI: ".oO[internalscanrange]Oo.",  #which might overwrite "CMS_zz4l_fai1"
    }

    repmap = {
              "cardstocombine": " ".join(["hzz4l_{}S_{}_{}.lumi{:.2f}.txt".format(channel, category, production.year, float(Luminosity(lumitype, production))) for channel, category, production in product(usechannels, usecategories, productions)] + alsocombine),
              "combinecardsfile": "hzz4l_4l.oO[combinecardsappend]Oo..txt",
              "workspacefile": "workspace.oO[workspacefileappend]Oo..root",
              "filename": "higgsCombine_.oO[append]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo...oO[method]Oo..mH125.root",
              "logfile": "log_.oO[method]Oo._.oO[append]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo..oO[expectfaiappend]Oo..txt",
              "expectedappend": "exp_.oO[expectfai]Oo.",
              "totallumi": "{:.2f}".format(totallumi),
              "observedappend": "obs",
              "setparameters":
                ",".join([
                  "CMS_zz4l_fai{}=.oO[expectfai]Oo.".format(analysis.fais.index(scanfai)+1)
                ] + [
                  "{}={}".format(k, v) for k, v in freeze.iteritems() if v is not None
                ]),
              "parameterranges": ":".join("{}={}".format(k, v) for k, v in parameterranges.iteritems()),
              "moreappend": moreappend,
              "turnoff": " ".join(turnoff),
              "workspacefileappend": workspacefileappend,
              "combinecardsappend": combinecardsappend,
              "algo": algo,
              "method": method,
              "robustfit": str(int(robustfit)),
              "setPOI": "" if POI==defaultPOI else "-P .oO[POI]Oo.",
              "POI": POI,
              "floatotherpois": str(int(floatothers if len(faiorder) > 1 else not fixfai)),
              "pointindex": "",
              "sqrts": ",".join("{:d}".format(_) for _ in sqrts),
              "physicsmodel": None,
              "physicsoptions": None,
              "morecombineoptions": morecombineoptions[method],
              "ntoys": str(ntoys),
              "saveorloadworkspace": "",
              "combine": "combineTool.py",
              "impactsstep1": "--doInitialFit",
              "impactsstep2": "--doFits",
              "impactsstep3": "-o .oO[filename]Oo.",
              "saveasdir": saveasdir,
              "fais": (
                " ".join("--PO {}".format(_) for _ in faiorder)
                + " "
                + " ".join("--PO {}asPOI{}".format(_, "" if str(_) == str(scanfai) else "relative") for _ in faiorder[:-1])
                + (" --PO {}asPOI".format(faiorder[-1]) if str(faiorder[-1]) != "fa1" else "")
              ),
              "freeze": ",".join(freeze),
              "freezeparameters": "--freezeParameters=.oO[freeze]Oo." if freeze else "",
              "setparametersforgrid": setparametersforgrid if setparametersforgrid else ".oO[setparameters]Oo.",
              "savespecifiednuis": "--saveSpecifiedNuis everything_but_binbybin" if usesystematics else "",
             }

    if method == "Impacts":
        repmap.update({
          "filename": ".oO[impactsfilename.oO[impactsstep]Oo.]Oo.",
          "logfile": repmap["logfile"].replace(".oO[method]Oo.", ".oO[method]Oo..oO[impactsstep]Oo."),
          "impactsfilename1": "higgsCombine_initialFit__.oO[append]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo..MultiDimFit.mH125.oO[seedforimpacts]Oo..root",
          "impactsfilename2": "higgsCombine_paramFit__.oO[append]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo._.oO[nuisanceforimpacts]Oo..MultiDimFit.mH125.oO[seedforimpacts]Oo..root",
          "impactsfilename3": "impacts_.oO[append]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo..json",
        })
        

    repmap["physicsmodel"] = "HiggsAnalysis.CombinedLimit.SpinZeroStructure:hzzAnomalousCouplingsFromHistograms"
    repmap["physicsoptions"] = "--PO sqrts=.oO[sqrts]Oo. --PO verbose --PO allowPMF .oO[fais]Oo. .oO[linearsystematics]Oo."
    repmap["savemu"] = "--saveSpecifiedFunc=" + ",".join(["CMS_zz4l_fa1"] + ["CMS_zz4l_fai"+str(i) for i, fai in enumerate(analysis.fais, start=1) if fai != scanfai] + ["fa3_ggH", "fCP_Htt", "RV", "RF"])
    repmap["setPOI"] = "-P CMS_zz4l_fai{}".format(analysis.fais.index(scanfai)+1)
    if analysis.isdecayonly:
        repmap["savemu"] = repmap["savemu"].replace(",fa3_ggH,fCP_Htt,RV", "")

    if analysis.isEFT:
        repmap["physicsmodel"] = "HiggsAnalysis.CombinedLimit.SpinZeroStructure:hzzAnomalousCouplingsFromHistogramsAi"
        repmap["savemu"] = "--saveSpecifiedFunc=g1,g2,g4,g1prime2,a2gg,a3gg,kappa,kappa_tilde"

    folder = os.path.join(config.repositorydir, "scans", subdirectory, "cards_{}".format(fullfoldername))
    mkdir_p(folder)
    with cd(folder):
        with open(".gitignore", "w") as f:
            f.write("*")
        datacards = makeDCsandWSs(productions, usecategories, usechannels, analysis, lumitype)
        for filename in alsocombine:
            for _ in filename, filename.replace(".input.root", ".txt"):
                if not os.path.exists(_): raise ValueError("{} does not exist!".format(_))
        tfiles.clear()
        with utilities.OneAtATime(replaceByMap(".oO[combinecardsfile]Oo..tmp", repmap), 5, task="running combineCards"):
            if not os.path.exists(replaceByMap(".oO[combinecardsfile]Oo.", repmap)):
                try:
                    subprocess.check_call(replaceByMap(combinecardstemplate, repmap), shell=True)
                except:
                    try:
                        raise
                    finally:
                        try:
                            os.remove(replaceByMap(".oO[combinecardsfile]Oo.", repmap))
                        except subprocess.CalledProcessError:
                            pass

        with open(replaceByMap(".oO[combinecardsfile]Oo.", repmap)) as f:
            for line in f:
                if line.startswith("shapesystematics group = "):
                    repmap["linearsystematics"] = "--PO linearsystematic:" + ",".join(line.split("=")[1].split())
                if re.match("CMS_scale_j *param", line):
                    repmap["parameterranges"] += ":CMS_scale_j=-1.0,1.0"

        with utilities.OneAtATime(replaceByMap(".oO[workspacefile]Oo..tmp", repmap), 5, task="running text2workspace"):
            if not os.path.exists(replaceByMap(".oO[workspacefile]Oo.", repmap)):
                subprocess.check_call(replaceByMap(createworkspacetemplate, repmap), shell=True)

        if onlyworkspace: return

        jobids = set()
        finalfiles = []

        if config.unblindscans and runobs:
          minimum = (float("nan"), float("inf"))
          plotlimitskwargs["combinelogs"].append([])
          for scanrange, internalscanrange in izip(scanranges, internalscanranges):
              repmap_obs = repmap.copy()
              repmap_obs.update({
                "npoints": str(scanrange[0]),
                "internalnpoints": str(internalscanrange[0]),
                "scanrange": "{},{}".format(*scanrange[1:3]),
                "internalscanrange": "{},{}".format(*internalscanrange[1:3]),
                "pointoffset": str(internalscanrange[3]),
                "scanrangeappend": ("" if scanrange==defaultscanrange else "_.oO[npoints]Oo.,.oO[scanrange]Oo.")+".oO[pointindex]Oo.",
                "expectfai": "0.0",  #starting point
                "append": ".oO[observedappend]Oo.",
                "expectfaiappend": "",
                "exporobs": "obs",
                "-t -1": "",
                "seedforimpacts": "",
              })

              jobids.add(runscan(repmap_obs, submitjobs=submitjobs))
              repmap_obs["impactsstep"] = "3"
              finalfiles.append(replaceByMap(".oO[filename]Oo.", repmap_obs))
              if not submitjobs and method == "MultiDimFit":
                  f = ROOT.TFile(replaceByMap(".oO[filename]Oo.", repmap_obs))
                  t = f.limit
                  for entry in t:
                      if t.deltaNLL+t.nll+t.nll0 < minimum[1]:
                          minimum = (getattr(t, POI), t.deltaNLL+t.nll+t.nll0)

              plotlimitskwargs["combinelogs"][-1].append(replaceByMap(".oO[logfile]Oo.", repmap_obs))

              if method == "FitDiagnostics":
                command = replaceByMap(diffnuisancescommand, repmap_obs)
                subprocess.check_call(command, shell=True)

              if method == "Impacts":
                command = replaceByMap(plotimpactscommand, repmap_obs)
                subprocess.check_call(command, shell=True)
                plotbasename = replaceByMap(".oO[saveasdir]Oo./impacts_.oO[append]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo.", repmap_obs)
                plotcopier.copy(plotbasename+".png")
                plotcopier.copy(plotbasename+".pdf")
                utilities.writeplotinfo(
                  plotbasename+".txt",
                  plotcopier=plotcopier,
                )

          minimum = minimum[0]
          del f

        for expectfai in expectvalues:
          plotlimitskwargs["combinelogs"].append([])
          for scanrange, internalscanrange in izip(scanranges, internalscanranges):
              repmap_exp = repmap.copy()
              if expectfai == "minimum":
                  expectfai = minimum
              repmap_exp.update({
                "npoints": str(scanrange[0]),
                "internalnpoints": str(internalscanrange[0]),
                "scanrange": "{},{}".format(*scanrange[1:]),
                "internalscanrange": "{},{}".format(*internalscanrange[1:3]),
                "pointoffset": str(internalscanrange[3]),
                "scanrangeappend": ("" if scanrange==defaultscanrange else "_.oO[npoints]Oo.,.oO[scanrange]Oo.")+".oO[pointindex]Oo.",
                "expectfai": str(expectfai),
                "append": ".oO[expectedappend]Oo.",
                "expectfaiappend": "_.oO[expectfai]Oo.",
                "exporobs": "exp",
                "-t -1": "-t .oO[ntoys]Oo.",
                "seedforimpacts": ""
              })
              jobids.add(runscan(repmap_exp, submitjobs=submitjobs))
              repmap_exp["impactsstep"] = "3"
              finalfiles.append(replaceByMap(".oO[filename]Oo.", repmap_exp))
              plotlimitskwargs["combinelogs"][-1].append(replaceByMap(".oO[logfile]Oo.", repmap_exp))

          if method == "FitDiagnostics":
            command = replaceByMap(diffnuisancescommand, repmap_exp)
            subprocess.check_call(command, shell=True)

          if method == "Impacts":
            command = replaceByMap(plotimpactscommand, repmap_exp)
            subprocess.check_call(command, shell=True)
            utilities.writeplotinfo(replaceByMap(".oO[saveasdir]Oo./impacts_.oO[append]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo..txt", repmap_exp))

        if None in jobids: jobids.remove(None)
        if jobids:
            submitjob("echo done", waitids=jobids, interactive=True, jobtime="0:0:10", jobname="wait")

        if ntry < maxntries:
            anydontexist = False
            for filename in finalfiles:
                if not utilities.existsandvalid(filename, "limit"): anydontexist = True
            if anydontexist:
                runcombine(*inputargs, **inputkwargs)
                return

        plotscans = []
        if config.unblindscans and runobs:
            plotscans.append("obs")
        for expectfai in expectvalues:
            if expectfai == "minimum":
                expectfai = minimum
            plotscans.append(expectfai)
        for ext in "png eps root pdf".split():
            plotname = plotname.replace("."+ext, "")
        plotname += replaceByMap(".oO[moreappend]Oo.", repmap)
        if scanranges != [defaultscanrange]:
            plotname += "".join("_{},{},{}".format(*scanrange) for scanrange in sorted(scanranges))

        if is2dscan:
            plotlimitsfunction = plotlimits2D
        else:
            plotlimitsfunction = plotlimits

        if faifor != "decay":
            plotname += "_"+faifor

        if method == "MultiDimFit" and scanranges:
            plotlimitsfunction(os.path.join(saveasdir, plotname), analysis, *plotscans, productions=productions, legendposition=legendposition, CLtextposition=CLtextposition, moreappend=replaceByMap(".oO[moreappend]Oo.", repmap), luminosity=totallumi, scanranges=scanranges, POI=POI, fixfai=fixfai, drawCMS=drawCMS, CMStext=CMStext, scanfai=scanfai, faifor=faifor, xaxislimits=xaxislimits, plotcopier=plotcopier, **plotlimitskwargs)
            for nuisance in plotnuisances:
                if nuisance == defaultPOI and fixfai: continue
                nuisanceplotname = plotname.replace("limit", plottitle(nuisance, analysis=analysis))
                plotlimitsfunction(os.path.join(saveasdir, nuisanceplotname), analysis, *plotscans, productions=productions, legendposition=legendposition, CLtextposition=CLtextposition, moreappend=replaceByMap(".oO[moreappend]Oo.", repmap), luminosity=totallumi, scanranges=scanranges, nuisance=nuisance, POI=POI, fixfai=fixfai, drawCMS=drawCMS, CMStext=CMStext, scanfai=scanfai, faiorder=faiorder, faifor=faifor, xaxislimits=xaxislimits, plotcopier=plotcopier, **plotlimitskwargs)

            utilities.writeplotinfo(
              os.path.join(saveasdir, plotname+".txt"),
              " ".join(pipes.quote(_) for _ in (
                "python", "limits.py", str(analysis), foldername,
                "--plotname="+plotname, "--poi="+POI, "--subdirectory="+subdirectory
              )),
              plotcopier=plotcopier,
            )

def runpermutations(analysis, foldername, *args, **kwargs):
  assert "faiorder" not in kwargs and "plotnuisances" not in kwargs
  scanfai = kwargs["scanfai"]

  analysis = Analysis(analysis)
  fullfoldername = "{}_{}".format(analysis, foldername)
  subdirectory = kwargs.get("subdirectory", "")

  for permutation in permutations([str(_) for _ in analysis.fais] + ["fa1"]):
    newkwargs = kwargs.copy()
    if permutation[0] != scanfai: continue

    newkwargs["faiorder"] = ",".join(permutation)

    plotnuisances = ["CMS_zz4l_fai"+str(i) for i in xrange(1, len(analysis.fais)+1)] + ["CMS_zz4l_fa1"]
    del plotnuisances[analysis.fais.index(scanfai)]
    newkwargs["plotnuisances"] = ",".join(plotnuisances)

    folder = os.path.join(config.repositorydir, "scans", subdirectory, "cards_{}".format(fullfoldername))
    mkdir_p(folder)

    with utilities.KeepWhileOpenFile(os.path.join(folder, "permutation_scan{}_{}.tmp".format(scanfai, newkwargs["faiorder"]))) as kwof:
      if not kwof: continue
      runcombine(analysis, foldername, *args, **newkwargs)

  newkwargs = kwargs.copy()
  newkwargs["floatothers"] = False
  with utilities.KeepWhileOpenFile(os.path.join(config.repositorydir, "scans", subdirectory, "cards_{}".format(fullfoldername), "permutation_scan{}_{}.tmp".format(scanfai, "fixothers"))) as kwof:
    if not kwof: return
    runcombine(analysis, foldername, *args, **newkwargs)

def rungrid(analysis, foldername, npoints, minpoint, maxpoint, *args, **kwargs):
  assert "setparametersforgrid" not in kwargs and "plotnuisances" not in kwargs
  analysis = Analysis(analysis)
  fullfoldername = "{}_{}".format(analysis, foldername)
  subdirectory = kwargs.get("subdirectory", "")
  scanfai = kwargs["scanfai"]

  otherfais = ["CMS_zz4l_fai"+str(i) for i in xrange(1, len(analysis.fais)+1)] + ["CMS_zz4l_fa1"]

  indicestoremove = [analysis.fais.index(scanfai)]
  if "faiorder" in kwargs:
    faiorder = tuple(Analysis(_) if _ != "fa1" else _ for _ in kwargs["faiorder"].split(","))
  else: 
    faiorder = tuple(analysis.fais) + ("fa1",)

  if str(faiorder[-1]) == "fa1":
    indicestoremove.append(len(otherfais)-1)
  else:
    assert False #not set up to handle minus sign
    indicestoremove.append(analysis.fais.index(faiorder[-1]))
  for _ in sorted(indicestoremove, reverse=True):
    del otherfais[_]

  otherfais = [_+"_relative" for _ in otherfais]

  if len(npoints) == len(minpoint) == len(maxpoint) == 1:
    npoints *= len(otherfais)
    minpoint *= len(otherfais)
    maxpoint *= len(otherfais)
  assert len(npoints) == len(minpoint) == len(maxpoint) == len(otherfais)

  plotnuisances = ["CMS_zz4l_fai"+str(i) for i in xrange(1, len(analysis.fais)+1)] + ["CMS_zz4l_fa1"]
  del plotnuisances[analysis.fais.index(scanfai)]
  kwargs["plotnuisances"] = ",".join(plotnuisances)

  for indices in product(*(xrange(npointsi+1) for npointsi in npoints)):
    otherfaivalues = [minpointi + (maxpointi-minpointi)*index/npointsi if npointsi != 0 else minpointi for index, npointsi, minpointi, maxpointi in izip_longest(indices, npoints, minpoint, maxpoint)]
    newkwargs = kwargs.copy()

    newkwargs["setparametersforgrid"] = ",".join("{}:{:g}".format(otherfai, otherfaivalue) for otherfai, otherfaivalue in izip(otherfais, otherfaivalues))

    folder = os.path.join(config.repositorydir, "scans", subdirectory, "cards_{}".format(fullfoldername))
    mkdir_p(folder)

    with utilities.KeepWhileOpenFile(os.path.join(folder, "grid_scan{}_{}.tmp".format(scanfai, newkwargs["setparametersforgrid"]))) as kwof:
      if not kwof: continue
      runcombine(analysis, foldername, *args, **newkwargs)

ntry = 0
maxntries = 3

def main():
    if sys.argv[1] == "runscan":
        function = runscan
        repmap = json.loads(sys.argv[2])
        args = [repmap, False]
        startkwargsfrom = 3
    elif sys.argv[1] == "runpermutations":
        function = runpermutations
        analysis = Analysis(sys.argv[2])
        foldername = sys.argv[3]
        args = [analysis, foldername]
        startkwargsfrom = 4
    elif sys.argv[1] == "rungrid":
        function = rungrid
        analysis = Analysis(sys.argv[2])
        foldername = sys.argv[3]
        npoints = [int(_) for _ in sys.argv[4].split(",")]
        minpoint = [float(_) for _ in sys.argv[5].split(",")]
        maxpoint = [float(_) for _ in sys.argv[6].split(",")]
        args = [analysis, foldername, npoints, minpoint, maxpoint]
        startkwargsfrom = 7
    else:
        function = runcombine
        analysis = Analysis(sys.argv[1])
        foldername = sys.argv[2]
        args = [analysis, foldername]
        startkwargsfrom = 3

    kwargs = {}
    for arg in sys.argv[startkwargsfrom:]:
        kw, kwarg = arg.split("=")
        if kw in kwargs:
            raise TypeError("Duplicate kwarg {}!".format(kw))
        kwargs[kw] = kwarg

    with PlotCopier() as pc:
        if function != runscan:
            assert "plotcopier" not in kwargs
            kwargs["plotcopier"] = pc
        function(*args, **kwargs)

if __name__ == "__main__":
    main()
