#!/usr/bin/env python
from Alignment.OfflineValidation.TkAlAllInOneTool.helperFunctions import replaceByMap  #easiest place to get it
from helperstuff import config, utilities
from helperstuff.combinehelpers import Luminosity
from helperstuff.datacard import makeDCsandWSs
from helperstuff.enums import Analysis, categories, Category, Channel, channels, Production, ProductionMode
from helperstuff.plotlimits import plotlimits, plottitle
from helperstuff.utilities import tfiles
from itertools import product
import os
import pipes
import ROOT
import shutil
import subprocess
import sys

combinecardstemplate = r"""
eval $(scram ru -sh) &&
combineCards.py .oO[cardstocombine]Oo. > .oO[combinecardsfile]Oo.
"""

createworkspacetemplate = r"""
eval $(scram ru -sh) &&
unbuffer text2workspace.py -m 125 .oO[combinecardsfile]Oo. -P HiggsAnalysis.CombinedLimit.SpinZeroStructure:multiSignalSpinZeroHiggs \
                           --PO verbose --PO allowPMF .oO[fixmu]Oo. -o .oO[workspacefile]Oo. -v 7 .oO[turnoff]Oo. \
                           |& tee log.text2workspace.oO[workspacefileappend]Oo. &&
exit ${PIPESTATUS[0]}
"""
runcombinetemplate = r"""
eval $(scram ru -sh) &&
combine -M MultiDimFit .oO[workspacefile]Oo. --algo=.oO[algo]Oo. --robustFit=.oO[robustfit]Oo. --points .oO[npoints]Oo. \
        --setPhysicsModelParameterRanges .oO[physicsmodelparameterranges]Oo. -m 125 .oO[setPOI]Oo. --floatOtherPOIs=1 \
        -n $1_.oO[append]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo. \
        .oO[-t -1]Oo. --setPhysicsModelParameters .oO[setphysicsmodelparameters]Oo. -V -v 3 --saveNLL \
        -S .oO[usesystematics]Oo. .oO[savemu]Oo. --saveSpecifiedNuis all --saveInactivePOI=1 |& tee log.oO[expectfaiappend]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo...oO[exporobs]Oo.
"""

def check_call_test(*args, **kwargs):
    print args[0]
#uncomment this for testing purposes
#subprocess.check_call = check_call_test

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
    expectmuffH = expectmuVVH = 1
    fixmuV = fixmuf = False
    plotnuisances = []
    defaultalgo = algo = "grid"
    robustfit = False
    defaultPOI = POI = "CMS_zz4l_fai1"
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
        elif kw == "expectmuffH":
            expectmuffH = float(kwarg)
        elif kw == "expectmuVVH":
            expectmuVVH = float(kwarg)
        elif kw == "fixmuV":
            fixmuV = bool(int(kwarg))
        elif kw == "fixmuf":
            fixmuf = bool(int(kwarg))
        elif kw == "plotnuisances":
            plotnuisances += kwarg.split(",")
        elif kw == "plotmus":
            if bool(int(kwarg)):
                plotnuisances += ["r_VVH", "r_ffH"]
        elif kw == "algo":
            algo = kwarg
        elif kw == "robustfit":
            robustfit = bool(int(kwarg))
        elif kw == "POI":
            POI = kwarg
        elif kw == "runobs":
            runobs = bool(int(kwarg))
        else:
            raise TypeError("Unknown kwarg: {}".format(kw))

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
        workspacefileappend += "_fixmuV"
    if fixmuf:
        workspacefileappend += "_fixmuf"
    if expectmuffH != 1:
        moreappend += "_muffH{}".format(expectmuffH)
    if expectmuVVH != 1:
        moreappend += "_muVVH{}".format(expectmuVVH)
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

    analysis = Analysis(analysis)
    foldername = "{}_{}".format(analysis, foldername)
    totallumi = sum(float(Luminosity(p, lumitype)) for p in productions)

    physicsmodelparameterranges = {
                                   "CMS_zz4l_fai1": "-1,1",
                                   POI: ".oO[scanrange]Oo.",  #which might overwrite "CMS_zz4l_fai1"
                                  }

    repmap = {
              "cardstocombine": " ".join("hzz4l_{}S_{}_{}.lumi{}.txt".format(channel, category, production.year, float(Luminosity(lumitype, production))) for channel, category, production in product(usechannels, usecategories, productions)),
              "combinecardsfile": "hzz4l_4l.oO[combinecardsappend]Oo..txt",
              "workspacefile": "workspace.oO[workspacefileappend]Oo..root",
              "filename": "higgsCombine_.oO[append]Oo..oO[moreappend]Oo..oO[scanrangeappend]Oo..MultiDimFit.mH125.root",
              "expectedappend": "exp_.oO[expectfai]Oo.",
              "totallumi": str(totallumi),
              "observedappend": "obs",
              "setphysicsmodelparameters": ".oO[expectmus,]Oo.CMS_zz4l_fai1=.oO[expectfai]Oo.",
              "physicsmodelparameterranges": ":".join("{}={}".format(k, v) for k, v in physicsmodelparameterranges.iteritems()),
              "usesystematics": str(int(usesystematics)),
              "moreappend": moreappend,
              "turnoff": " ".join(turnoff),
              "workspacefileappend": workspacefileappend,
              "combinecardsappend": combinecardsappend,
              "expectmus,": ("" if fixmuf else "r_ffH=.oO[expectmuffH]Oo.,") + ("" if fixmuV else "r_VVH=.oO[expectmuVVH]Oo.,"),
              "expectmuffH": str(expectmuffH),
              "expectmuVVH": str(expectmuVVH),
              "fixmu": ("--PO muVFixed" if fixmuV else "") + ("--PO mufFixed" if fixmuf else ""),
              "savemu": "--saveSpecifiedFunc=" + ",".join(mu for mu, fix in (("r_VVH", fixmuV), ("r_ffH", fixmuf)) if not fix and mu!=POI),
              "algo": algo,
              "robustfit": str(int(robustfit)),
              "setPOI": "" if POI==defaultPOI else "-P .oO[POI]Oo.",
              "POI": POI,
             }
    folder = os.path.join(config.repositorydir, "CMSSW_7_6_5/src/HiggsAnalysis/HZZ4l_Combination/CreateDatacards", subdirectory, "cards_{}".format(foldername))
    utilities.mkdir_p(folder)
    with utilities.cd(folder):
        #must make it for all categories and channels even if not using them all because of mu definition!
        makeDCsandWSs(productions, categories, channels, analysis, lumitype)
        tfiles.clear()
        with open(".gitignore", "w") as f:
            f.write("*")
        with utilities.OneAtATime(replaceByMap(".oO[combinecardsfile]Oo..tmp", repmap), 5, task="running combineCards"):
            if not os.path.exists(replaceByMap(".oO[combinecardsfile]Oo.", repmap)):
                subprocess.check_call(replaceByMap(combinecardstemplate, repmap), shell=True)
        with utilities.OneAtATime(replaceByMap(".oO[workspacefile]Oo..tmp", repmap), 5, task="running text2workspace"):
            if not os.path.exists(replaceByMap(".oO[workspacefile]Oo.", repmap)):
                subprocess.check_call(replaceByMap(createworkspacetemplate, repmap), shell=True)

        if config.unblindscans and runobs:
          minimum = (float("nan"), float("inf"))
          for scanrange in scanranges:
              repmap_obs = repmap.copy()
              repmap_obs.update({
                "npoints": str(scanrange[0]),
                "scanrange": "{},{}".format(*scanrange[1:]),
                "scanrangeappend": "" if scanrange==defaultscanrange else "_.oO[npoints]Oo.,.oO[scanrange]Oo.",
                "expectfai": "0.0",  #starting point
                "append": ".oO[observedappend]Oo.",
                "expectfaiappend": "",
                "exporobs": "obs",
                "-t -1": "",
              })
              if not os.path.exists(replaceByMap(".oO[filename]Oo.", repmap_obs)):
                with utilities.OneAtATime(replaceByMap(".oO[filename]Oo.", repmap_obs)+".tmp", 30):
                  subprocess.check_call(replaceByMap(runcombinetemplate, repmap_obs), shell=True)
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
                "scanrangeappend": "" if scanrange==defaultscanrange else "_.oO[npoints]Oo.,.oO[scanrange]Oo.",
                "expectfai": str(expectfai),
                "append": ".oO[expectedappend]Oo.",
                "expectfaiappend": "_.oO[expectfai]Oo.",
                "exporobs": "exp",
                "-t -1": "-t -1",
              })
              if not os.path.exists(replaceByMap(".oO[filename]Oo.", repmap_exp)):
                with utilities.OneAtATime(replaceByMap(".oO[filename]Oo.", repmap_exp)+".tmp", 30):
                  subprocess.check_call(replaceByMap(runcombinetemplate, repmap_exp), shell=True)

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
            if nuisance!=POI:
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
            

if __name__ == "__main__":
    analysis = Analysis(sys.argv[1])
    foldername = sys.argv[2]
    kwargs = {}
    for arg in sys.argv[3:]:
        kw, kwarg = arg.split("=")
        if kw in kwargs:
            raise TypeError("Duplicate kwarg {}!".format(kw))
        kwargs[kw] = kwarg
    runcombine(analysis, foldername, **kwargs)
