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
combine -M MultiDimFit .oO[workspacefile]Oo. --algo=grid --points 100 -m 125 -n $1_exp -t -1 --setPhysicsModelParameters r=1,CMS_zz4l_fg4=0.0 --expectSignal=1 -V -v 3 |& tee log.exp
"""
observationcombinetemplate = """
eval $(scram ru -sh) &&
combine -M MultiDimFit .oO[workspacefile]Oo. --algo=grid --points 100 -m 125 -n $1_obs       --setPhysicsModelParameters r=1,CMS_zz4l_fg4=0.0                  -V -v 3 |& tee log.obs
"""


def runcombine(analysis, foldername, **kwargs):
    usechannels = channels
    observation = False
    for kw, kwarg in kwargs.iteritems():
        if kw == "channels":
            usechannels = [Channel(c) for c in kwarg.split(",")]
        elif kw == "observation":
            observation = bool(int(kwarg))
        else:
            raise TypeError("Unknown kwarg: {}".format(kw))


    analysis = Analysis(analysis)
    foldername = "{}_{}".format(analysis, foldername)
    repmap = {
              "foldername": pipes.quote(foldername),
              "analysis": str(analysis),
              "cardstocombine": " ".join("hzz4l_{}S_8TeV.txt".format(channel) for channel in usechannels),
              "workspacefile": "floatMu.root",
              "expectedfilename": "higgsCombine_exp.MultiDimFit.mH125.root",
              "observedfilename": "higgsCombine_obs.MultiDimFit.mH125.root",
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
                            if observation:
                                rates = getrates(channel, "fordata", config.productionforcombine)
                            else:
                                rates = getrates(channel, "forexpectedscan")
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
            if not os.path.exists(repmap["expectedfilename"]):
                subprocess.check_call(replaceByMap(runcombinetemplate, repmap), shell=True)
            if observation and not os.path.exists(repmap["observedfilename"]):
                subprocess.check_call(replaceByMap(observationcombinetemplate, repmap), shell=True)
            saveasdir = os.path.join(config.plotsbasedir, "limits", foldername)
            try:
                os.makedirs(saveasdir)
            except OSError:
                pass
            kwargs = {}
            if observation:
                kwargs.update(observedfilename = repmap["observedfilename"], production = config.productionforcombine)
            plotlimits(repmap["expectedfilename"], "CMS_zz4l_fg4", os.path.join(saveasdir, "limit"), analysis.title(), **kwargs)

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
