#!/usr/bin/env python

from Alignment.OfflineValidation.TkAlAllInOneTool.helperFunctions import replaceByMap  #easiest place to get it
from helperstuff import config
from helperstuff import filemanager
from helperstuff.combinehelpers import getrates, Luminosity
from helperstuff.enums import Analysis, Channel, channels, Production
from helperstuff.plotlimits import plotlimits
from helperstuff.replacesystematics import replacesystematics
from itertools import product
import os
import pipes
import ROOT
import subprocess
import sys

makeworkspacestemplate = """
eval $(scram ru -sh) &&
python make_prop_DCsandWSs.py -i SM_inputs_8TeV -a .oO[foldername]Oo. -A .oO[analysis]Oo. -P .oO[production]Oo. --sigmaVVai=T1:.oO[SM_XS_VV]Oo.,T2:.oO[BSM_XS_VV]Oo.,T4:.oO[int_XS_VV]Oo. --newMu
"""

createworkspacetemplate = """
eval $(scram ru -sh) &&
combineCards.py .oO[cardstocombine]Oo. > hzz4l_4l.txt &&
text2workspace.py -m 125 hzz4l_4l.txt -P HiggsAnalysis.CombinedLimit.SpinZeroStructure:multiSignalSpinZeroHiggs --PO allowPMF --PO sqrts=13 .oO[morePO]Oo. -o .oO[workspacefile]Oo. -v 7 |& tee log.text2workspace &&
exit ${PIPESTATUS[0]}
"""
runcombinetemplate = """
eval $(scram ru -sh) &&
combine -M MultiDimFit .oO[workspacefile]Oo. --algo=grid --points .oO[npoints]Oo. --setPhysicsModelParameterRanges CMS_zz4l_fai1=.oO[scanrange]Oo. -m 125 -n $1_.oO[append]Oo..oO[moreappend]Oo. -t -1 --setPhysicsModelParameters r=1,CMS_zz4l_fai1=.oO[expectfai]Oo. --expectSignal=1 -V -v 3 --saveNLL -S .oO[usesystematics]Oo. |& tee log_.oO[expectfai]Oo..oO[moreappend]Oo..exp
"""
observationcombinetemplate = """
eval $(scram ru -sh) &&
combine -M MultiDimFit .oO[workspacefile]Oo. --algo=grid --points .oO[npoints]Oo. --setPhysicsModelParameterRanges CMS_zz4l_fai1=.oO[scanrange]Oo. -m 125 -n $1_.oO[append]Oo..oO[moreappend]Oo.       --setPhysicsModelParameters r=1,CMS_zz4l_fai1=.oO[expectfai]Oo.                  -V -v 3 --saveNLL -S .oO[usesystematics]Oo. |& tee log.oO[moreappend]Oo..obs
"""

def check_call_test(*args, **kwargs):
    print args[0]
#uncomment this for testing purposes
#subprocess.check_call = check_call_test

def runcombine(analysis, foldername, **kwargs):
    usechannels = channels
    expectvalues = [0.0]
    plotname = "limit"
    legendposition = (.2, .7, .6, .9)
    CLtextposition = "left"
    productions = config.productionsforcombine
    usesystematics = True
    usebkg = True
    runobs = True
    defaultscanrange = scanrange = (-1.0, 1.0)
    defaultnpoints = npoints = 100
    defaultscenario = scenario = 1
    scalemuvmuftogether = False

    for kw, kwarg in kwargs.iteritems():
        if kw == "channels":
            usechannels = [Channel(c) for c in kwarg.split(",")]
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
        elif kw == "usebkg":
            usebkg = bool(int(kwarg))
        elif kw == "scanrange":
            try:
                scanrange = tuple(float(a) for a in kwarg.split(","))
                if len(scanrange) != 2: raise ValueError
            except ValueError:
                raise ValueError("scanrange has to contain 2 floats separated by a comma!")
        elif kw == "npoints":
            npoints = int(kwarg)
        elif kw == "scenario":
            scenario = int(kwarg)
        elif kw == "scalemuvmuftogether":
            scalemuvmuftogether = bool(int(kwarg))
        elif kw == "runobs":
            runobs = bool(int(kwarg))
        else:
            raise TypeError("Unknown kwarg: {}".format(kw))

    if not config.unblindscans:
        luminosity = float(Luminosity("forexpectedscan"))
    else:
        luminosity = None

    years = [p.year for p in productions]
    if len(set(years)) != len(years):
        raise ValueError("Some of your productions are from the same year!")

    workspacefileappend = ""
    if scalemuvmuftogether:
        workspacefileappend += "_scalemuvmuftogether"
    print workspacefileappend

    moreappend = ".oO[workspacefileappend]Oo."
    if not usesystematics:
        moreappend += "_nosystematics"
    if not (npoints == defaultnpoints and scanrange == defaultscanrange):
        moreappend += "_{},{},{}".format(npoints, *scanrange)

    analysis = Analysis(analysis)
    foldername = "{}_{}".format(analysis, foldername)
    if scenario != defaultscenario:
        foldername += "_scenario{}".format(scenario)
    repmap = {
              "foldername": pipes.quote(foldername),
              "analysis": str(analysis),
              "cardstocombine": " ".join("hzz4l_{}S_{}.txt".format(channel, production.year) for channel, production in product(usechannels, productions)),
              "workspacefile": "floatMu.oO[workspacefileappend]Oo..root",
              "filename": "higgsCombine_.oO[append]Oo..oO[moreappend]Oo..MultiDimFit.mH125.root",
              "expectedappend": "exp_.oO[expectfai]Oo.",
              "observedappend": "obs",
              "workspacefileappend": workspacefileappend,
              "usesystematics": str(int(usesystematics)),
              "moreappend": moreappend,
              "npoints": str(npoints),
              "scanrange": "{},{}".format(*scanrange),
              "SM_XS_VV": str(analysis.SM_XS_VV),
              "BSM_XS_VV": str(analysis.BSM_XS_VV),
              "int_XS_VV": str(analysis.int_XS_VV),
              "morePO": "--PO scalemuvmuftogether" if scalemuvmuftogether else "",
             }
    with filemanager.cd(os.path.join(config.repositorydir, "CMSSW_7_6_5/src/HiggsAnalysis/HZZ4l_Combination/CreateDatacards")):
        for production in productions:
            production = Production(production)
            if not all(os.path.exists("cards_{}/HCG/125/{}/hzz4l_{}S_{}.input.root".format(foldername, production.year, channel, production.year)) for channel in usechannels):
                makeworkspacesmap = repmap.copy()
                makeworkspacesmap["production"] = str(production)
                subprocess.check_call(replaceByMap(makeworkspacestemplate, makeworkspacesmap), shell=True)
        with open("cards_{}/.gitignore".format(foldername), "w") as f:
            f.write("*")
        with filemanager.cd("cards_{}/HCG/125/{}".format(foldername, production.year)):
            #replace rates
            for channel, production in product(channels, productions):
                if channel in usechannels:
                    with open("hzz4l_{}S_{}.txt".format(channel, production.year)) as f:
                        contents = f.read()
                        if "\n#rate" in contents: continue #already did this
                    for line in contents.split("\n"):
                        if line.startswith("rate"):
                            if config.unblindscans:
                                rates = getrates(channel, "fordata", production, usebkg=usebkg)
                            else:
                                rates = getrates(channel, "forexpectedscan", production, usebkg=usebkg)
                            contents = contents.replace(line, "#"+line+"\n"+rates)
                            break
                    with open("hzz4l_{}S_{}.txt".format(channel, production.year), "w") as f:
                        f.write(contents)
                else:
                    if os.path.exists("hzz4l_{}S_{}.txt".format(channel, production.year)):
                        os.remove("hzz4l_{}S_{}.txt".format(channel, production.year))
                        os.remove("hzz4l_{}S_{}.input.root".format(channel, production.year))
            if not os.path.exists(repmap["workspacefile"]):
                for channel, production in product(usechannels, productions):
                    replacesystematics(channel, production, scenario=scenario, luminosity=luminosity)
                subprocess.check_call(replaceByMap(createworkspacetemplate, repmap), shell=True)

            if config.unblindscans and runobs:
                repmap_obs = repmap.copy()
                repmap_obs["expectfai"] = "0.0"  #starting point
                repmap_obs["append"] = ".oO[observedappend]Oo."
                if not os.path.exists(replaceByMap(".oO[filename]Oo.", repmap_obs)):
                    subprocess.check_call(replaceByMap(observationcombinetemplate, repmap_obs), shell=True)
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
                if not os.path.exists(replaceByMap(".oO[filename]Oo.", repmap_exp)):
                    subprocess.check_call(replaceByMap(runcombinetemplate, repmap_exp), shell=True)

            saveasdir = os.path.join(config.plotsbasedir, "limits", foldername)
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
            plotname += moreappend
            plotname = replaceByMap(plotname, repmap)
            plotlimits(os.path.join(saveasdir, plotname), analysis, *plotscans, productions=productions, legendposition=legendposition, CLtextposition=CLtextposition, moreappend=replaceByMap(moreappend, repmap))
            with open(os.path.join(saveasdir, plotname+".txt"), "w") as f:
                f.write(" ".join(["python"]+sys.argv))
                f.write("\n\n\n")
                f.write("python limits.py ")
                for arg in sys.argv[1:]:
                    if "=" in arg: continue
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
