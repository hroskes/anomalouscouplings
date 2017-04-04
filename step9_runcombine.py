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
from helperstuff.datacard import makeDCsandWSs
from helperstuff.enums import Analysis, categories, Category, Channel, channels, Production, ProductionMode
from helperstuff.plotlimits import plotlimits, plottitle
from helperstuff.submitjob import submitjob
from helperstuff.utilities import tfiles

combinecardstemplate = r"""
eval $(scram ru -sh) &&
combineCards.py .oO[cardstocombine]Oo. > .oO[combinecardsfile]Oo.
"""

createworkspacetemplate = r"""
eval $(scram ru -sh) &&
unbuffer text2workspace.py -m 125 .oO[combinecardsfile]Oo. -P HiggsAnalysis.CombinedLimit.SpinZeroStructure:multiSignalSpinZeroHiggs \
                           --PO sqrts=.oO[sqrts]Oo. --PO verbose --PO allowPMF -o .oO[workspacefile]Oo. -v 7 .oO[turnoff]Oo. \
                           |& tee log.text2workspace.oO[workspacefileappend]Oo. &&
exit ${PIPESTATUS[0]}
"""
runcombinetemplate = r"""
eval $(scram ru -sh) &&
combine -M MultiDimFit .oO[workspacefile]Oo. --algo=.oO[algo]Oo. --robustFit=.oO[robustfit]Oo. --points .oO[npoints]Oo. \
        --setPhysicsModelParameterRanges .oO[physicsmodelparameterranges]Oo. -m 125 .oO[setPOI]Oo. --floatOtherPOIs=.oO[floatotherpois]Oo. \
        -n $1_.oO[append]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo. .oO[selectpoints]Oo. \
        --includePOIEdges=1 --X-rtd OPTIMIZE_BOUNDS=0 --X-rtd TMCSO_AdaptivePseudoAsimov=0 \
        .oO[-t -1]Oo. --setPhysicsModelParameters .oO[setphysicsmodelparameters]Oo. -V -v 3 --saveNLL \
        -S .oO[usesystematics]Oo. .oO[savemu]Oo. --saveSpecifiedNuis all --saveInactivePOI=1 |& tee log.oO[expectfaiappend]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo...oO[exporobs]Oo.
"""

def check_call_test(*args, **kwargs):
    print args[0]
#uncomment this for testing purposes
#subprocess.check_call = check_call_test

def runscan(repmap, submitjobs, directory=None):
  if submitjobs:
    npoints = int(repmap["npoints"])
    jobids = set()
    individualfilenames = set()

    repmap_final = repmap.copy()
    repmap_final["selectpoints"] = ""
    finalfilename = replaceByMap(".oO[filename]Oo.", repmap_final)
    if os.path.exists(finalfilename): return

    for i in range(npoints+1):
      repmap_i = repmap.copy()
      repmap_i.update({
        "selectpoints": "--firstPoint .oO[pointindex]Oo. --lastPoint .oO[pointindex]Oo.",
        "pointindex": str(i),
      })
      replaceByMap(runcombinetemplate, repmap_i) #sanity check that replaceByMap works
      individualfilenames.add(replaceByMap(".oO[filename]Oo.", repmap_i))
      if os.path.exists(replaceByMap(".oO[filename]Oo.", repmap_i)):
        continue
      job = "{} runscan {} directory={}".format(
                                                os.path.join(config.repositorydir, "./step9_runcombine.py"),
                                                pipes.quote(json.dumps(repmap_i)),
                                                pipes.quote(os.getcwd()),
                                               )
      jobname = replaceByMap(".oO[expectfaiappend]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo...oO[exporobs]Oo.", repmap_i).lstrip("_")
      jobid = submitjob(job, jobname=jobname, jobtime="1-0:0:0", morerepmap=repmap_i)
      jobids.add(jobid)

    haddjobid = submitjob(" ".join(["hadd", finalfilename] + list(individualfilenames)), jobname="hadd", jobtime="1:0:0", waitids=jobids, docd=True)
    return haddjobid

  else:
    cwd = os.path.realpath(os.getcwd())
    if directory is None: directory = cwd
    directory = os.path.realpath(directory)
    if cwd != directory:
#      if utilities.LSB_JOBID() is None:
#        raise ValueError("Should call runscan from directory except in a batch job")
      shutil.copytree(directory, os.path.basename(directory.rstrip("/")), symlinks=True)
      cdto = os.path.basename(directory.rstrip("/"))
    else:
      cdto = "."

    repmap = repmap.copy()
    if "selectpoints" not in repmap:
      repmap["selectpoints"] = ""
    filename = replaceByMap(".oO[filename]Oo.", repmap)
    tmpfile = os.path.join(directory, filename+".tmp")
    if not os.path.exists(filename):
           #utilities.KeepWhileOpenFile(tmpfile, message=utilities.LSB_JOBID()), \
      with utilities.cd(cdto), \
           utilities.LSF_creating(os.path.join(directory, filename)):
        subprocess.check_call(replaceByMap(runcombinetemplate, repmap), shell=True)


def runcombine(analysis, foldername, **kwargs):
    usechannels = channels
    usecategories = categories
    expectvalues = [0.0]
    plotname = "limit"
    legendposition = (.2, .7, .6, .9)
    CLtextposition = "left"
    productions = config.productionsforcombine
    usesystematics = True
    runobs = True
    subdirectory = ""
    defaultscanrange = (100, -1.0, 1.0)
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
    if config.unblindscans:
        lumitype = "fordata"
    else:
        lumitype = "forexpectedscan"
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
            if config.unblindscans:
                raise TypeError("For unblindscans, if you want to adjust the luminosity do it in the Production class (in enums.py)")
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
        else:
            raise TypeError("Unknown kwarg: {}".format(kw))

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
        moreappend += "_scan"+plottitle(POI)
    if fixfai:
        moreappend += "_fixfai"
    if alsocombine:
        combinecardsappend += "_" + alsocombinename
        if sqrts is None:
            raise ValueError("Have to provide sqrts if you provide alsocombine!")

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
              "cardstocombine": " ".join(["hzz4l_{}S_{}_{}.lumi{}.txt".format(channel, category, production.year, float(Luminosity(lumitype, production))) for channel, category, production in product(usechannels, usecategories, productions)] + alsocombine),
              "combinecardsfile": "hzz4l_4l.oO[combinecardsappend]Oo..txt",
              "workspacefile": "workspace.oO[workspacefileappend]Oo..root",
              "filename": "higgsCombine_.oO[append]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo..MultiDimFit.mH125.root",
              "expectedappend": "exp_.oO[expectfai]Oo.",
              "totallumi": str(totallumi),
              "observedappend": "obs",
              "setphysicsmodelparameters": "CMS_zz4l_fai1=.oO[expectfai]Oo.",
              "physicsmodelparameterranges": ":".join("{}={}".format(k, v) for k, v in physicsmodelparameterranges.iteritems()),
              "usesystematics": str(int(usesystematics)),
              "moreappend": moreappend,
              "turnoff": " ".join(turnoff),
              "workspacefileappend": workspacefileappend,
              "combinecardsappend": combinecardsappend,
              "savemu": "--saveSpecifiedFunc=" + ",".join(mu for mu, fix in (("muV,muV_scaled", fixmuV), ("muf,muf_scaled", fixmuf)) if not fix and mu!=POI),
              "algo": algo,
              "robustfit": str(int(robustfit)),
              "setPOI": "" if POI==defaultPOI else "-P .oO[POI]Oo.",
              "POI": POI,
              "floatotherpois": str(int(not fixfai)),
              "pointindex": "",
              "sqrts": ",".join("{:d}".format(_) for _ in sqrts),
             }
    folder = os.path.join(config.repositorydir, "CMSSW_7_6_5/src/HiggsAnalysis/HZZ4l_Combination/CreateDatacards", subdirectory, "cards_{}".format(foldername))
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
                subprocess.check_call(replaceByMap(combinecardstemplate, repmap), shell=True)
        with utilities.OneAtATime(replaceByMap(".oO[workspacefile]Oo..tmp", repmap), 5, task="running text2workspace"):
            if not os.path.exists(replaceByMap(".oO[workspacefile]Oo.", repmap)):
                subprocess.check_call(replaceByMap(createworkspacetemplate, repmap), shell=True)

        jobids = set()

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

        if None in jobids: jobids.remove(None)
        if jobids:
            submitjob("echo done", waitids=jobids, interactive=True, jobtime="0:0:10", jobname="wait")

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
        plotlimits(os.path.join(saveasdir, plotname), analysis, *plotscans, productions=productions, legendposition=legendposition, CLtextposition=CLtextposition, moreappend=replaceByMap(".oO[moreappend]Oo.", repmap), luminosity=totallumi, scanranges=scanranges, POI=POI)
        for nuisance in plotnuisances:
            if nuisance == POI: continue
            if nuisance == "CMS_zz4l_fai1" and fixfai: continue
            nuisanceplotname = plotname.replace("limit", plottitle(nuisance))
            plotlimits(os.path.join(saveasdir, nuisanceplotname), analysis, *plotscans, productions=productions, legendposition=legendposition, CLtextposition=CLtextposition, moreappend=replaceByMap(".oO[moreappend]Oo.", repmap), luminosity=totallumi, scanranges=scanranges, nuisance=nuisance, POI=POI)

    with open(os.path.join(saveasdir, plotname+".txt"), "w") as f:
        f.write(" ".join(["python"]+[pipes.quote(_) for _ in sys.argv]))
        f.write("\n\n\n")
        f.write("python limits.py ")
        for arg in sys.argv[1:]:
            if "=" in arg and "subdirectory=" not in arg: continue
            f.write(pipes.quote(arg)+" ")
        f.write("plotname="+plotname+" ")
        f.write("\n\n\n\n\n\ngit info:\n\n")
        f.write(subprocess.check_output(["git", "rev-parse", "HEAD"]))
        f.write("\n")
        f.write(subprocess.check_output(["git", "status"]))
        f.write("\n")
        f.write(subprocess.check_output(["git", "diff"]))
            

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
