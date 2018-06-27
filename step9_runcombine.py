#!/usr/bin/env python
import glob
from itertools import product
import json
import os
import pipes
import shutil
import subprocess
import sys

import ROOT

from Alignment.OfflineValidation.TkAlAllInOneTool.helperFunctions import replaceByMap  #easiest place to get it

from helperstuff import config, utilities
from helperstuff.combinehelpers import Luminosity
from helperstuff.copyplots import copyplots
from helperstuff.datacard import makeDCsandWSs
from helperstuff.enums import Analysis, categories, Category, Channel, channels, Production, ProductionMode
from helperstuff.plotlimits import plotlimits, plotlimits2D, plottitle
from helperstuff.submitjob import submitjob
from helperstuff.utilities import requirecmsenv, tfiles

requirecmsenv(os.path.join(config.repositorydir, "CMSSW_8_1_0"))

combinecardstemplate = r"""
combineCards.py .oO[cardstocombine]Oo. > .oO[combinecardsfile]Oo.
"""

createworkspacetemplate = r"""
unbuffer text2workspace.py -m 125 .oO[combinecardsfile]Oo. -P .oO[physicsmodel]Oo. \
                           .oO[physicsoptions]Oo. -o .oO[workspacefile]Oo. -v 7 .oO[turnoff]Oo. \
                           |& tee log.text2workspace.oO[workspacefileappend]Oo. &&
exit ${PIPESTATUS[0]}
"""
runcombinetemplate = r"""
combine -M MultiDimFit .oO[workspacefile]Oo. --algo=.oO[algo]Oo. --robustFit=.oO[robustfit]Oo. --points .oO[npoints]Oo. \
        --setParameterRanges .oO[physicsmodelparameterranges]Oo. -m 125 .oO[setPOI]Oo. --floatOtherPOIs=.oO[floatotherpois]Oo. \
        -n $1_.oO[append]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo. .oO[selectpoints]Oo. \
        --alignEdges=1 --X-rtd OPTIMIZE_BOUNDS=0 --X-rtd TMCSO_AdaptivePseudoAsimov=0 \
        .oO[-t -1]Oo. --setParameters .oO[setphysicsmodelparameters]Oo. -V -v 3 --saveNLL \
        -S .oO[usesystematics]Oo. .oO[savemu]Oo. --saveSpecifiedNuis all --saveInactivePOI=1 |& tee log.oO[expectfaiappend]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo...oO[exporobs]Oo.
"""

def check_call_test(*args, **kwargs):
    print args[0]
#uncomment this for testing purposes
#subprocess.check_call = subprocess.check_output = check_call_test

def runscan(repmap, submitjobs, directory=None):
  if submitjobs:
    utilities.mkdir_p("jobs")
    npoints = int(repmap["npoints"])
    jobids = set()
    individualfilenames = set()

    repmap_final = repmap.copy()
    repmap_final["selectpoints"] = ""
    finalfilename = replaceByMap(".oO[filename]Oo.", repmap_final)
    if os.path.exists(finalfilename): return

    for i in range(npoints):
      repmap_i = repmap.copy()
      repmap_i.update({
        "selectpoints": "--firstPoint .oO[pointindex]Oo. --lastPoint .oO[pointindex]Oo.",
        "pointindex": str(i),
      })
      replaceByMap(runcombinetemplate, repmap_i) #sanity check that replaceByMap works
      individualfilenames.add(replaceByMap("jobs/.oO[filename]Oo.", repmap_i))

      ############################
      #compatibility
      if os.path.exists(replaceByMap(".oO[filename]Oo.", repmap_i)):
        shutil.move(replaceByMap(".oO[filename]Oo.", repmap_i), "jobs/")
      ############################

      if os.path.exists(replaceByMap("jobs/.oO[filename]Oo.", repmap_i)):
        try:
          f = ROOT.TFile(replaceByMap("jobs/.oO[filename]Oo.", repmap_i))
          f.limit
        except Exception as e:
          try:
            del f
            os.remove(replaceByMap("jobs/.oO[filename]Oo.", repmap_i))
          except:
            pass
        else:
          del f

      if (os.path.exists(replaceByMap("jobs/.oO[filename]Oo.", repmap_i))
       or os.path.exists(replaceByMap("jobs/.oO[filename]Oo..tmp", repmap_i))):
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
        "morerepmap": repmap_i
      }
      if config.host == "MARCC":
        submitjobkwargs["queue"] = "lrgmem"
        submitjobkwargs["memory"] = "10G"
      jobid = submitjob(job, **submitjobkwargs)
      jobids.add(jobid)

    haddjobid = submitjob(" ".join([os.path.join(config.repositorydir, "helperstuff", "hadd.py"), finalfilename] + list(individualfilenames)), jobname="hadd", jobtime="1:0:0", waitids=jobids, docd=True)
    return haddjobid

  else:
    cwd = os.path.realpath(os.getcwd())
    if directory is None: directory = cwd
    directory = os.path.realpath(directory)
    if cwd != directory:
#      if utilities.LSB_JOBID() is None:
#        raise ValueError("Should call runscan from directory except in a batch job")
      shutil.copy(os.path.join(directory, replaceByMap(".oO[workspacefile]Oo.", repmap)), ".")
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
    filename = replaceByMap(".oO[filename]Oo.", repmap)
    tmpfile = os.path.join(directory, filename+".tmp")
    logfile = replaceByMap("log.oO[expectfaiappend]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo...oO[exporobs]Oo.", repmap)
    if not os.path.exists(filename):
      with utilities.cd(cdto), \
           utilities.OneAtATime(tmpfile, 30, message=utilities.LSB_JOBID()), \
           utilities.LSF_creating(os.path.join(directory, filename), skipifexists=True), \
           utilities.LSF_creating(os.path.join(directory, logfile), skipifexists=True):
        if not os.path.exists(filename):
          subprocess.check_call(replaceByMap(runcombinetemplate, repmap), shell=True)


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
    defaultscanrange = (101, -1.0, 1.0)
    scanranges = [defaultscanrange]
    defaultusesignalproductionmodes = usesignalproductionmodes = {ProductionMode(p) for p in ("ggH", "VBF", "ZH", "WH", "ttH")}
    usebkg = True
    fixmuV = fixmuf = fixfai = False
    plotnuisances = []
    defaultalgo = algo = "grid"
    robustfit = False
    defaultPOI = POI = "CMS_zz4l_fai1"
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
            if analysis.dimensions != 2:
                raise ValueError("scanfai is only for 2D analyses")
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
        else:
            raise TypeError("Unknown kwarg: {}".format(kw))

    if runobs and not config.unblindscans:
        raise TypeError("Can't unblind scans!")
    if runobs and lumitype != "fordata":
        raise TypeError("For unblindscans, if you want to adjust the luminosity do it in the Production class (in enums.py)")
    if scanfai != analysis and scanfai not in analysis.fais:
        raise ValueError("scanfai for {} has to be ".format(analysis) + " or ".join(str(_) for _ in list(analysis.fais)+[analysis]))

    is2dscan = (scanfai == analysis and analysis.dimensions == 2)
    if is2dscan:
      scanranges = list((points**2, min, max) for points, min, max in scanranges)

    if submitjobs:
        if not utilities.inscreen():
            raise RuntimeError("submitjobs should be run from a screen session!")
        if "minimum" in expectvalues:
            raise ValueError("Can't run scan for minimum using submitjobs")

    if len(productions) > 1 and lumitype != "fordata":
        raise TypeError("If there's >1 production, have to use lumitype == 'fordata'")

    years = [p.year for p in productions]
    if len(set(years)) != len(years):
        raise ValueError("Some of your productions are from the same year!")

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
    if set(usechannels) != set(channels):
        combinecardsappend += "_" + ",".join(sorted(str(c) for c in usechannels))
    if set(usecategories) != set(categories):
        combinecardsappend += "_" + ",".join(sorted(str(c) for c in usecategories))
    if not usebkg:
        workspacefileappend += "_nobkg"
        turnoff.append("--PO nobkg")
    if not usesystematics:
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
    if alsocombine:
        combinecardsappend += "_" + alsocombinename
        if sqrts is None:
            raise ValueError("Have to provide sqrts if you provide alsocombine!")
    if analysis.dimensions == 2 and not is2dscan:
        workspacefileappend += "_scan{}".format(scanfai)

    if set(usecategories) != {Category("Untagged")} and analysis.isdecayonly:
        raise ValueError("For decay only analysis have to specify categories=Untagged")
    if set(usechannels) != {Channel("2e2mu")} and config.LHE:
        raise ValueError("For LHE analysis have to specify channels=2e2mu")

    if sqrts is None:
        sqrts = [13]

    analysis = Analysis(analysis)
    foldername = "{}_{}".format(analysis, foldername)
    totallumi = sum(float(Luminosity(p, lumitype)) for p in productions)

    physicsmodelparameterranges = {
                                   "CMS_zz4l_fai1": "-1,1",
                                   POI: ".oO[scanrange]Oo.",  #which might overwrite "CMS_zz4l_fai1"
                                  }

    repmap = {
              "cardstocombine": " ".join(["hzz4l_{}S_{}_{}.lumi{:.2f}.txt".format(channel, category, production.year, float(Luminosity(lumitype, production))) for channel, category, production in product(usechannels, usecategories, productions)] + alsocombine),
              "combinecardsfile": "hzz4l_4l.oO[combinecardsappend]Oo..txt",
              "workspacefile": "workspace.oO[workspacefileappend]Oo..root",
              "filename": "higgsCombine_.oO[append]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo..MultiDimFit.mH125.root",
              "expectedappend": "exp_.oO[expectfai]Oo.",
              "totallumi": "{:.2f}".format(totallumi),
              "observedappend": "obs",
              "setphysicsmodelparameters": "CMS_zz4l_fai1=.oO[expectfai]Oo.",
              "physicsmodelparameterranges": ":".join("{}={}".format(k, v) for k, v in physicsmodelparameterranges.iteritems()),
              "usesystematics": str(int(usesystematics)),
              "moreappend": moreappend,
              "turnoff": " ".join(turnoff),
              "workspacefileappend": workspacefileappend,
              "combinecardsappend": combinecardsappend,
              "algo": algo,
              "robustfit": str(int(robustfit)),
              "setPOI": "" if POI==defaultPOI else "-P .oO[POI]Oo.",
              "POI": POI,
              "floatotherpois": str(int(not fixfai)),
              "pointindex": "",
              "sqrts": ",".join("{:d}".format(_) for _ in sqrts),
              "physicsmodel": None,
              "physicsoptions": None,
             }
    if analysis.usehistogramsforcombine:
        repmap["physicsmodel"] = "HiggsAnalysis.CombinedLimit.SpinZeroStructure:hzzAnomalousCouplingsFromHistograms"
        repmap["physicsoptions"] = "--PO sqrts=.oO[sqrts]Oo. --PO verbose --PO allowPMF --PO .oO[analysis]Oo."
        repmap["savemu"] = "--saveSpecifiedFunc=" + ",".join(mu for mu, fix in (("muV", fixmuV), ("muf", fixmuf)) if not fix and mu!=POI and not analysis.isdecayonly)
    else:
        if analysis.dimensions == 2 and analysis.isdecayonly:
            repmap["physicsmodel"] = "HiggsAnalysis.CombinedLimit.SpinZeroStructure:spinZeroHiggs"
            repmap["physicsoptions"] = "--PO allowPMF"
            if scanfai == analysis: repmap["physicsoptions"] += " --PO fai2asPOI"
            elif scanfai == analysis.fais[0]: pass
            elif scanfai == analysis.fais[1]: repmap["physicsoptions"] += " --PO fai2asPOI --PO fai1fixed"
            repmap["savemu"] = ""
        elif analysis.dimensions == 2 and not analysis.isdecayonly:
            assert False
        elif analysis.dimensions != 2 and analysis.isdecayonly:
            assert False
        else:
            repmap["physicsmodel"] = "HiggsAnalysis.CombinedLimit.SpinZeroStructure:multiSignalSpinZeroHiggs"
            repmap["physicsoptions"] = "--PO sqrts=.oO[sqrts]Oo. --PO verbose --PO allowPMF"
            repmap["savemu"] = "--saveSpecifiedFunc=" + ",".join(mu for mu, fix in (("muV", fixmuV), ("muf", fixmuf)) if not fix and mu!=POI and not analysis.isdecayonly)

    folder = os.path.join(config.repositorydir, "scans", subdirectory, "cards_{}".format(foldername))
    utilities.mkdir_p(folder)
    with utilities.cd(folder):
        with open(".gitignore", "w") as f:
            f.write("*")
        #must make it for all categories and channels even if not using them all because of mu definition!
        makeDCsandWSs(productions, categories, channels, analysis, lumitype)
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
        with utilities.OneAtATime(replaceByMap(".oO[workspacefile]Oo..tmp", repmap), 5, task="running text2workspace"):
            if not os.path.exists(replaceByMap(".oO[workspacefile]Oo.", repmap)):
                subprocess.check_call(replaceByMap(createworkspacetemplate, repmap), shell=True)

        jobids = set()
        finalfiles = []

        if config.unblindscans and runobs:
          minimum = (float("nan"), float("inf"))
          for scanrange in scanranges:
              repmap_obs = repmap.copy()
              repmap_obs.update({
                "npoints": str(scanrange[0]),
                "scanrange": "{},{}".format(*scanrange[1:]),
                "scanrangeappend": ("" if scanrange==defaultscanrange else "_.oO[npoints]Oo.,.oO[scanrange]Oo.")+".oO[pointindex]Oo.",
                "expectfai": "0.0",  #starting point
                "append": ".oO[observedappend]Oo.",
                "expectfaiappend": "",
                "exporobs": "obs",
                "-t -1": "",
              })
              jobids.add(runscan(repmap_obs, submitjobs=submitjobs))
              finalfiles.append(replaceByMap(".oO[filename]Oo.", repmap_obs))
              if not submitjobs:
                  f = ROOT.TFile(replaceByMap(".oO[filename]Oo.", repmap_obs))
                  t = f.limit
                  for entry in t:
                      if t.deltaNLL+t.nll+t.nll0 < minimum[1]:
                          minimum = (getattr(t, POI), t.deltaNLL+t.nll+t.nll0)
          minimum = minimum[0]
          del f

        for scanrange in scanranges:
          for expectfai in expectvalues:
              repmap_exp = repmap.copy()
              if expectfai == "minimum":
                  expectfai = minimum
              repmap_exp.update({
                "npoints": str(scanrange[0]),
                "scanrange": "{},{}".format(*scanrange[1:]),
                "scanrangeappend": ("" if scanrange==defaultscanrange else "_.oO[npoints]Oo.,.oO[scanrange]Oo.")+".oO[pointindex]Oo.",
                "expectfai": str(expectfai),
                "append": ".oO[expectedappend]Oo.",
                "expectfaiappend": "_.oO[expectfai]Oo.",
                "exporobs": "exp",
                "-t -1": "-t -1",
              })
              jobids.add(runscan(repmap_exp, submitjobs=submitjobs))
              finalfiles.append(replaceByMap(".oO[filename]Oo.", repmap_exp))

        if None in jobids: jobids.remove(None)
        if jobids:
            submitjob("echo done", waitids=jobids, interactive=True, jobtime="0:0:10", jobname="wait")

        if ntry < maxntries:
            e = None
            for filename in finalfiles:
                try:
                    f = ROOT.TFile(filename)
                    f.limit
                except Exception as e:
                    try:
                        del f
                        os.remove(filename)
                    except:
                        pass
                else:
                    del f
            if e is not None:
                runcombine(*inputargs, **inputkwargs)
                return

        saveasdir = os.path.join(config.plotsbasedir, "limits", subdirectory, foldername)
        try:
            os.makedirs(saveasdir)
        except OSError:
            pass
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
        plotlimitsfunction(os.path.join(saveasdir, plotname), analysis, *plotscans, productions=productions, legendposition=legendposition, CLtextposition=CLtextposition, moreappend=replaceByMap(".oO[moreappend]Oo.", repmap), luminosity=totallumi, scanranges=scanranges, POI=POI, fixfai=fixfai, drawCMS=drawCMS, CMStext=CMStext, scanfai=scanfai, faifor=faifor, xaxislimits=xaxislimits)
        for nuisance in plotnuisances:
            if plottitle(nuisance) == plottitle(POI): continue
            if nuisance == "CMS_zz4l_fai1" and fixfai: continue
            nuisanceplotname = plotname.replace("limit", plottitle(nuisance))
            plotlimitsfunction(os.path.join(saveasdir, nuisanceplotname), analysis, *plotscans, productions=productions, legendposition=legendposition, CLtextposition=CLtextposition, moreappend=replaceByMap(".oO[moreappend]Oo.", repmap), luminosity=totallumi, scanranges=scanranges, nuisance=nuisance, POI=POI, fixfai=fixfai, drawCMS=drawCMS, CMStext=CMStext, faifor=faifor, xaxislimits=xaxislimits)

    with open(os.path.join(saveasdir, plotname+".txt"), "w") as f:
        f.write(" ".join(["python"]+[pipes.quote(_) for _ in sys.argv]))
        f.write("\n\n\n")
        f.write("python limits.py ")
        for arg in sys.argv[1:]:
            if "=" in arg and "subdirectory=" not in arg: continue
            f.write(pipes.quote(arg)+" ")
        f.write("--plotname="+plotname+" ")
        f.write("\n\n\n\n\n\ngit info:\n\n")
        f.write(subprocess.check_output(["git", "rev-parse", "HEAD"]))
        f.write("\n")
        f.write(subprocess.check_output(["git", "status"]))
        f.write("\n")
        f.write(subprocess.check_output(["git", "diff"]))

    copyplots(os.path.join("limits", subdirectory, foldername))

ntry = 0
maxntries = 3

def main():
    if sys.argv[1] == "runscan":
        function = runscan
        repmap = json.loads(sys.argv[2])
        args = [repmap, False]
    else:
        function = runcombine
        analysis = Analysis(sys.argv[1])
        foldername = sys.argv[2]
        args = [analysis, foldername]

    kwargs = {}
    for arg in sys.argv[3:]:
        kw, kwarg = arg.split("=")
        if kw in kwargs:
            raise TypeError("Duplicate kwarg {}!".format(kw))
        kwargs[kw] = kwarg

    function(*args, **kwargs)

if __name__ == "__main__":
    main()
