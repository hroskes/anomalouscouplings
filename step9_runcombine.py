from Alignment.OfflineValidation.TkAlAllInOneTool.helperFunctions import replaceByMap  #easiest place to get it
from helperstuff import config
from helperstuff import filemanager
from helperstuff.combinehelpers import getrates
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
python make_prop_DCsandWSs.py -i SM_inputs_8TeV -a .oO[foldername]Oo. -A .oO[analysis]Oo. -P .oO[production]Oo.
"""

createworkspacetemplate = """
eval $(scram ru -sh) &&
combineCards.py .oO[cardstocombine]Oo. > hzz4l_4l.txt &&
text2workspace.py -m 125 hzz4l_4l.txt -P HiggsAnalysis.CombinedLimit.SpinZeroStructure:spinZeroHiggs --PO=muFloating -o .oO[workspacefile]Oo. -v 7 |& tee log.text2workspace
"""
runcombinetemplate = """
eval $(scram ru -sh) &&
combine -M MultiDimFit .oO[workspacefile]Oo. --algo=grid --points 100 -m 125 -n $1_.oO[append]Oo. -t -1 --setPhysicsModelParameters r=1,CMS_zz4l_fg4=.oO[expectfai]Oo. --expectSignal=1 -V -v 3 --saveNLL |& tee log_.oO[expectfai]Oo..exp
"""
observationcombinetemplate = """
eval $(scram ru -sh) &&
combine -M MultiDimFit .oO[workspacefile]Oo. --algo=grid --points 100 -m 125 -n $1_.oO[append]Oo.       --setPhysicsModelParameters r=1,CMS_zz4l_fg4=.oO[expectfai]Oo.                  -V -v 3 --saveNLL |& tee log.obs
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
    for kw, kwarg in kwargs.iteritems():
        if kw == "channels":
            usechannels = [Channel(c) for c in kwarg.split(",")]
        elif kw == "expectvalues":
            try:
                expectvalues = [fai if fai=="minimum" else float(fai) for fai in kwarg.split(",")]
            except ValueError:
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
        else:
            raise TypeError("Unknown kwarg: {}".format(kw))

    years = [p.year for p in productions]
    if len(set(years)) != len(years):
        raise ValueError("Some of your productions are from the same year!")

    analysis = Analysis(analysis)
    foldername = "{}_{}".format(analysis, foldername)
    repmap = {
              "foldername": pipes.quote(foldername),
              "analysis": str(analysis),
              "cardstocombine": " ".join("hzz4l_{}S_{}.txt".format(channel, production.year) for channel, production in product(usechannels, productions)),
              "workspacefile": "floatMu.root",
              "filename": "higgsCombine_.oO[append]Oo..MultiDimFit.mH125.root",
              "expectedappend": "exp_.oO[expectfai]Oo.",
              "observedappend": "obs",
             }
    with filemanager.cd(os.path.join(config.repositorydir, "CMSSW_7_6_5/src/HiggsAnalysis/HZZ4l_Combination/CreateDatacards")):
        for production in productions:
            production = Production(production)
            if not all(os.path.exists("cards_{}/HCG/125/hzz4l_{}S_{}.input.root".format(foldername, channel, production.year)) for channel in channels):
                makeworkspacesmap = repmap.copy()
                makeworkspacesmap["production"] = str(production)
                subprocess.check_call(replaceByMap(makeworkspacestemplate, makeworkspacesmap), shell=True)
        with open("cards_{}/.gitignore".format(foldername), "w") as f:
            f.write("*")
        with filemanager.cd("cards_{}/HCG/125".format(foldername)):
            #replace rates
            for channel, production in product(channels, productions):
                if channel in usechannels:
                    with open("hzz4l_{}S_{}.txt".format(channel, production.year)) as f:
                        contents = f.read()
                        if "\n#rate" in contents: continue #already did this
                    for line in contents.split("\n"):
                        if line.startswith("rate"):
                            if config.unblindscans:
                                rates = getrates(channel, "fordata", production)
                            else:
                                rates = getrates(channel, "forexpectedscan", production)
                            contents = contents.replace(line, "#"+line+"\n"+rates)
                            break
                    with open("hzz4l_{}S_{}.txt".format(channel, production.year), "w") as f:
                        f.write(contents)
                else:
                    os.remove("hzz4l_{}S_{}.txt".format(channel, production.year))
                    os.remove("hzz4l_{}S_{}.input.root".format(channel, production.year))
            if not os.path.exists(repmap["workspacefile"]):
                for channel, production in product(usechannels, productions):
                    replacesystematics(channel, production)
                subprocess.check_call(replaceByMap(createworkspacetemplate, repmap), shell=True)

            if config.unblindscans:
                repmap_obs = repmap.copy()
                repmap_obs["expectfai"] = "0.0"  #starting point
                repmap_obs["append"] = ".oO[observedappend]Oo."
                if not os.path.exists(replaceByMap(".oO[filename]Oo.", repmap_obs)):
                    subprocess.check_call(replaceByMap(observationcombinetemplate, repmap_obs), shell=True)
                f = ROOT.TFile(replaceByMap(".oO[filename]Oo.", repmap_obs))
                f.limit.GetEntry(0)
                minimum = f.limit.CMS_zz4l_fg4
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
            if config.unblindscans:
                plotscans.append("obs")
            for expectfai in expectvalues:
                if expectfai == "minimum":
                    expectfai = minimum
                plotscans.append(expectfai)
            for ext in "png eps root pdf".split():
                plotname = plotname.replace("."+ext, "")
            plotlimits(os.path.join(saveasdir, plotname), analysis, *plotscans, productions=productions, legendposition=legendposition, CLtextposition=CLtextposition)
            with open(os.path.join(saveasdir, plotname+".txt"), "w") as f:
                f.write(" ".join(["python"]+sys.argv) + "\n")

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
