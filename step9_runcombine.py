from Alignment.OfflineValidation.TkAlAllInOneTool.helperFunctions import replaceByMap  #easiest place to get it
from helperstuff import config
from helperstuff import filemanager
from helperstuff.combinehelpers import getrates
from helperstuff.enums import Analysis, Channel, channels
from helperstuff.plotlimits import plotlimits
from helperstuff.replacesystematics import replacesystematics
import os
import pipes
import subprocess
import sys

makeworkspacestemplate = """
eval $(scram ru -sh) &&
python make_prop_DCsandWSs.py -i SM_inputs_8TeV -a .oO[foldername]Oo. -A .oO[analysis]Oo.
"""

createworkspacetemplate = """
eval $(scram ru -sh) &&
combineCards.py .oO[cardstocombine]Oo. > hzz4l_4l_8TeV.txt &&
text2workspace.py -m 125 hzz4l_4l_8TeV.txt -P HiggsAnalysis.CombinedLimit.SpinZeroStructure:spinZeroHiggs --PO=muFloating -o .oO[workspacefile]Oo. -v 7 |& tee log.text2workspace
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
    for kw, kwarg in kwargs.iteritems():
        if kw == "channels":
            usechannels = [Channel(c) for c in kwarg.split(",")]
        elif kw == "expectvalues":
            try:
                expectvalues = [float(fai) for fai in kwarg.split(",")]
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
        else:
            raise TypeError("Unknown kwarg: {}".format(kw))


    analysis = Analysis(analysis)
    foldername = "{}_{}".format(analysis, foldername)
    repmap = {
              "foldername": pipes.quote(foldername),
              "analysis": str(analysis),
              "cardstocombine": " ".join("hzz4l_{}S_8TeV.txt".format(channel) for channel in usechannels),
              "workspacefile": "floatMu.root",
              "filename": "higgsCombine_.oO[append]Oo..MultiDimFit.mH125.root",
              "expectedappend": "exp_.oO[expectfai]Oo.",
              "observedappend": "obs",
             }
    with filemanager.cd(os.path.join(config.repositorydir, "CMSSW_7_6_5/src/HiggsAnalysis/HZZ4l_Combination/CreateDatacards")):
        if not os.path.exists("cards_{}".format(foldername)):
            subprocess.check_call(replaceByMap(makeworkspacestemplate, repmap), shell=True)
        with open("cards_{}/.gitignore".format(foldername), "w") as f:
            f.write("*")
        with filemanager.cd("cards_{}/HCG/125".format(foldername)):
            #replace rates
            for channel in channels:
                if channel in usechannels:
                    with open("hzz4l_{}S_8TeV.txt".format(channel)) as f:
                        contents = f.read()
                        if "\n#rate" in contents: continue #already did this
                    for line in contents.split("\n"):
                        if line.startswith("rate"):
                            if config.unblindscans:
                                rates = getrates(channel, "fordata", config.productionforcombine)
                            else:
                                rates = getrates(channel, "forexpectedscan", config.productionforcombine)
                            contents = contents.replace(line, "#"+line+"\n"+rates)
                            break
                    with open("hzz4l_{}S_8TeV.txt".format(channel), "w") as f:
                        f.write(contents)
                else:
                    os.remove("hzz4l_{}S_8TeV.txt".format(channel))
                    os.remove("hzz4l_{}S_8TeV.input.root".format(channel))
            if not os.path.exists(repmap["workspacefile"]):
                for channel in usechannels:
                    replacesystematics(channel)
                subprocess.check_call(replaceByMap(createworkspacetemplate, repmap), shell=True)

            for expectfai in expectvalues:
                repmap_exp = repmap.copy()
                repmap_exp["expectfai"] = str(expectfai)
                repmap_exp["append"] = ".oO[expectedappend]Oo."
                if not os.path.exists(replaceByMap(".oO[filename]Oo.", repmap_exp)):
                    subprocess.check_call(replaceByMap(runcombinetemplate, repmap_exp), shell=True)

            if config.unblindscans:
                repmap_obs = repmap.copy()
                repmap_obs["expectfai"] = "0.0"  #starting point
                repmap_obs["append"] = ".oO[observedappend]Oo."
                if not os.path.exists(replaceByMap(".oO[filename]Oo.", repmap_obs)):
                    subprocess.check_call(replaceByMap(observationcombinetemplate, repmap_obs), shell=True)

            saveasdir = os.path.join(config.plotsbasedir, "limits", foldername)
            try:
                os.makedirs(saveasdir)
            except OSError:
                pass
            plotscans = []
            if config.unblindscans:
                plotscans.append("obs")
            plotscans += expectvalues
            for ext in "png eps root pdf".split():
                plotname = plotname.replace("."+ext, "")
            plotlimits(os.path.join(saveasdir, plotname), analysis, *plotscans, production=config.productionforcombine, legendposition=legendposition, CLtextposition=CLtextposition)
            with open(os.path.join(saveasdir, plotname+".txt"), "w") as f:
                f.write(" ".join(["python"]+sys.argv))

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
