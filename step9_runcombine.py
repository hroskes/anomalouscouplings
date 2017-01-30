#!/usr/bin/env python
from Alignment.OfflineValidation.TkAlAllInOneTool.helperFunctions import replaceByMap  #easiest place to get it
from helperstuff import config, utilities
from helperstuff.combinehelpers import Luminosity
from helperstuff.datacard import makeDCandWS
from helperstuff.enums import Analysis, categories, Category, Channel, channels, Production, ProductionMode
from helperstuff.plotlimits import plotlimits
from helperstuff.replacesystematics import replacesystematics
from itertools import product
import os
import pipes
import ROOT
import shutil
import subprocess
import sys

createworkspacetemplate = r"""
eval $(scram ru -sh) &&
combineCards.py .oO[cardstocombine]Oo. > .oO[combinecardsfile]Oo. &&
unbuffer text2workspace.py -m 125 .oO[combinecardsfile]Oo. -P HiggsAnalysis.CombinedLimit.SpinZeroStructure:multiSignalSpinZeroHiggs \
                           --PO verbose --PO=muFloating,allowPMF -o .oO[workspacefile]Oo. -v 7 .oO[turnoff]Oo. \
                           |& tee log.text2workspace.oO[workspacefileappend]Oo. &&
exit ${PIPESTATUS[0]}
"""
runcombinetemplate = r"""
eval $(scram ru -sh) &&
combine -M MultiDimFit .oO[workspacefile]Oo. --algo=grid --points .oO[npoints]Oo. \
        --setPhysicsModelParameterRanges CMS_zz4l_fai1=.oO[scanrange]Oo. -m 125 -n $1_.oO[append]Oo..oO[moreappend]Oo. \
        .oO[-t -1]Oo. --setPhysicsModelParameters .oO[setphysicsmodelparameters]Oo. -V -v 3 --saveNLL \
        -S .oO[usesystematics]Oo. --saveSpecifiedFunc=r_VVH,r_ggH |& tee log.oO[expectfaiappend]Oo..oO[moreappend]Oo...oO[exporobs]Oo.
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
    subdirectory = ""
    defaultscanrange = scanrange = (-1.0, 1.0)
    defaultnpoints = npoints = 100
    defaultusesignalproductionmodes = usesignalproductionmodes = (ProductionMode(p) for p in ("ggH", "VBF", "ZH", "WH"))
    usebkg = True
    expectmuggH = expectmuVVH = 1
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
        elif kw == "scanrange":
            try:
                scanrange = tuple(float(a) for a in kwarg.split(","))
                if len(scanrange) != 2: raise ValueError
            except ValueError:
                raise ValueError("scanrange has to contain 2 floats separated by a comma!")
        elif kw == "npoints":
            npoints = int(kwarg)
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
        elif kw == "expectmuggH":
            expectmuggH = float(kwarg)
        elif kw == "expectmuVVH":
            expectmuVVH = float(kwarg)
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
    if expectmuggH != 1:
        moreappend += "_muggH{}".format(expectmuggH)
    if expectmuVVH != 1:
        moreappend += "_muVVH{}".format(expectmuVVH)
    if usesignalproductionmodes != defaultusesignalproductionmodes:
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
    if not (npoints == defaultnpoints and scanrange == defaultscanrange):
        moreappend += "_{},{},{}".format(npoints, *scanrange)

    analysis = Analysis(analysis)
    foldername = "{}_{}".format(analysis, foldername)
    totallumi = sum(float(Luminosity(p, lumitype)) for p in productions)

    repmap = {
              "cardstocombine": " ".join("hzz4l_{}S_{}_{}.lumi{}.txt".format(channel, category, production.year, float(Luminosity(lumitype, production))) for channel, category, production in product(usechannels, usecategories, productions)),
              "combinecardsfile": "hzz4l_4l.oO[combinecardsappend]Oo..txt",
              "workspacefile": "floatMu.oO[workspacefileappend]Oo..root",
              "filename": "higgsCombine_.oO[append]Oo..oO[moreappend]Oo..MultiDimFit.mH125.root",
              "expectedappend": "exp_.oO[expectfai]Oo.",
              "totallumi": str(totallumi),
              "observedappend": "obs",
              "setphysicsmodelparameters": ".oO[expectrs]Oo.,CMS_zz4l_fai1=.oO[expectfai]Oo.",
              "usesystematics": str(int(usesystematics)),
              "moreappend": moreappend,
              "npoints": str(npoints),
              "scanrange": "{},{}".format(*scanrange),
              "turnoff": " ".join(turnoff),
              "workspacefileappend": workspacefileappend,
              "combinecardsappend": combinecardsappend,
              "expectrs": "r_ggH=.oO[expectmuggH]Oo.,r_VVH=.oO[expectmuVVH]Oo.",
              "expectmuggH": str(expectmuggH),
              "expectmuVVH": str(expectmuVVH),
             }
    folder = os.path.join(config.repositorydir, "CMSSW_7_6_5/src/HiggsAnalysis/HZZ4l_Combination/CreateDatacards", subdirectory, "cards_{}".format(foldername))
    utilities.mkdir_p(folder)
    with utilities.cd(folder):
        for production, category, channel in product(productions, usecategories, usechannels):
            makeDCandWS(production, category, channel, analysis, lumitype)
        with open(".gitignore", "w") as f:
            f.write("*")
        if not os.path.exists(replaceByMap(".oO[workspacefile]Oo.", repmap)):
            subprocess.check_call(replaceByMap(createworkspacetemplate, repmap), shell=True)

        if config.unblindscans:
            repmap_obs = repmap.copy()
            repmap_obs["expectfai"] = "0.0"  #starting point
            repmap_obs["append"] = ".oO[observedappend]Oo."
            repmap_obs["expectfaiappend"] = ""
            repmap_obs["exporobs"] = "obs"
            repmap_obs["-t -1"] = ""
            if not os.path.exists(replaceByMap(".oO[filename]Oo.", repmap_obs)):
                subprocess.check_call(replaceByMap(runcombinetemplate, repmap_obs), shell=True)
            f = ROOT.TFile(replaceByMap(".oO[filename]Oo.", repmap_obs))
            minimum = (float("nan"), float("inf"))
            for entry in f.limit:
                if f.limit.deltaNLL < minimum[1]:
                    minimum = (f.limit.CMS_zz4l_fai1, f.limit.deltaNLL)
            minimum = minimum[0]
            del f

        for expectfai in expectvalues:
            repmap_exp = repmap.copy()
            if expectfai == "minimum":
                expectfai = minimum
            repmap_exp["expectfai"] = str(expectfai)
            repmap_exp["append"] = ".oO[expectedappend]Oo."
            repmap_exp["expectfaiappend"] = "_.oO[expectfai]Oo."
            repmap_exp["exporobs"] = "exp"
            repmap_exp["-t -1"] = "-t -1"
            if not os.path.exists(replaceByMap(".oO[filename]Oo.", repmap_exp)):
                subprocess.check_call(replaceByMap(runcombinetemplate, repmap_exp), shell=True)

        saveasdir = os.path.join(config.plotsbasedir, "limits", subdirectory, foldername)
        try:
            os.makedirs(saveasdir)
        except OSError:
            pass
        plotscans = []
        if config.unblindscans:
            plotscans.append("obs")
        for expectfai in expectvalues:
            if expectfai == "minimum":
                expectfai = minimum
            plotscans.append(expectfai)
        for ext in "png eps root pdf".split():
            plotname = plotname.replace("."+ext, "")
        plotname += replaceByMap(".oO[moreappend]Oo.", repmap)
        plotlimits(os.path.join(saveasdir, plotname), analysis, *plotscans, productions=productions, legendposition=legendposition, CLtextposition=CLtextposition, moreappend=replaceByMap(".oO[moreappend]Oo.", repmap), luminosity=totallumi)
        with open(os.path.join(saveasdir, plotname+".txt"), "w") as f:
            f.write(" ".join(["python"]+sys.argv))
            f.write("\n\n\n")
            f.write("python limits.py ")
            for arg in sys.argv[1:]:
                if "=" in arg and "subdirectory=" not in arg: continue
                f.write(arg+" ")
            f.write("plotname="+plotname+" ")
            f.write("\n")

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
