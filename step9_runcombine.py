from Alignment.OfflineValidation.TkAlAllInOneTool.helperFunctions import replaceByMap  #easiest place to get it
from helperstuff import config
from helperstuff import filemanager
from helperstuff.combinehelpers import getrates
from helperstuff.enums import Analysis, Channel, channels
from helperstuff.plotlimits import plotlimits
import os
import pipes
import subprocess
import sys

makeworkspacestemplate = """
eval $(scram ru -sh) &&
python make_prop_DCsandWSs.py -i SM_inputs_8TeV -a .oO[foldername]Oo. -A .oO[analysis]Oo.
"""

runcombinetemplate = """
eval $(scram ru -sh) &&
combineCards.py .oO[cardstocombine]Oo. > hzz4l_4l_8TeV.txt &&
text2workspace.py -m 125 hzz4l_4l_8TeV.txt -P HiggsAnalysis.CombinedLimit.SpinZeroStructure:spinZeroHiggs --PO=muFloating -o fixedMu_8TeV.root -v 7 &&
combine -M MultiDimFit fixedMu_8TeV.root --algo=grid --points 100 -m 125 -n $1_8TeV -t -1 --setPhysicsModelParameters CMS_zz4l_fg4=0.0 --expectSignal=1 -V -v 3
"""

def runcombine(analysis, foldername, **kwargs):
    usechannels = channels
    for kw, kwarg in kwargs.iteritems():
        if kw == "channels":
            usechannels = [Channel(c) for c in kwarg.split(",")]
        else:
            raise ValueError("Unknown kwarg: {}".format(kw))


    analysis = Analysis(analysis)
    foldername = "{}_{}".format(analysis, foldername)
    repmap = {
              "foldername": pipes.quote(foldername),
              "analysis": str(analysis),
              "cardstocombine": " ".join("hzz4l_{}S_8TeV.txt".format(channel) for channel in usechannels)
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
                    for line in contents.split("\n"):
                        if line.startswith("rate"):
                            contents = contents.replace(line, "#"+line+"\n"+getrates(channel))
                            break
                    with open("hzz4l_{}S_8TeV.txt".format(channel), "w") as f:
                        f.write(contents)
                else:
                    os.remove("hzz4l_{}S_8TeV.txt".format(channel))
                    os.remove("hzz4l_{}S_8TeV.input.root".format(channel))
            if not os.path.exists("higgsCombine_8TeV.MultiDimFit.mH125.root"):
                subprocess.check_call(replaceByMap(runcombinetemplate, repmap), shell=True)
            saveasdir = os.path.join(config.plotsbasedir, "limits", foldername)
            try:
                os.makedirs(saveasdir)
            except OSError:
                pass
            plotlimits("higgsCombine_8TeV.MultiDimFit.mH125.root", "CMS_zz4l_fg4", os.path.join(saveasdir, "limit"), analysis.title())

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
