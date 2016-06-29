from Alignment.OfflineValidation.TkAlAllInOneTool.helperFunctions import replaceByMap  #easiest place to get it
from helperstuff import config
from helperstuff import filemanager
from helperstuff.combinehelpers import getrates
from helperstuff.enums import channels
from helperstuff.plotlimits import plotlimits
import os
import pipes
import subprocess
import sys

makeworkspacestemplate = """
eval $(scram ru -sh) &&
python make_prop_DCsandWSs.py -i SM_inputs_8TeV -a .oO[foldername]Oo.
"""

runcombinetemplate = """
eval $(scram ru -sh) &&
combineCards.py hzz4l_2e2muS_8TeV.txt hzz4l_4muS_8TeV.txt hzz4l_4eS_8TeV.txt > hzz4l_4l_8TeV.txt &&
text2workspace.py -m 125 hzz4l_4l_8TeV.txt -P HiggsAnalysis.CombinedLimit.SpinZeroStructure:spinZeroHiggs --PO=muFloating -o fixedMu_8TeV.root -v 7 &&
combine -M MultiDimFit fixedMu_8TeV.root --algo=grid --points 100 -m 125 -n $1_8TeV -t -1 --setPhysicsModelParameters CMS_zz4l_fg4=0.0 --expectSignal=1 -V -v 3
"""

def runcombine(foldername):
    repmap = {
              "foldername": pipes.quote(foldername),
             }
    with filemanager.cd(os.path.join(config.repositorydir, "CMSSW_7_6_5/src/HiggsAnalysis/HZZ4l_Combination/CreateDatacards")):
        subprocess.check_call(replaceByMap(makeworkspacestemplate, repmap), shell=True)
        with open("cards_{}/.gitignore".format(foldername), "w") as f:
            f.write("*")
        with filemanager.cd("cards_{}/HCG/125".format(foldername)):
            #replace rates
            for channel in channels:
                with open("hzz4l_{}S_8TeV.txt".format(channel)) as f:
                    contents = f.read()
                for line in contents.split("\n"):
                    if line.startswith("rate"):
                        contents = contents.replace(line, "#"+line+"\n"+getrates(channel))
                        break
                with open("hzz4l_{}S_8TeV.txt".format(channel), "w") as f:
                     f.write(contents)
            subprocess.check_call(replaceByMap(runcombinetemplate, repmap), shell=True)
            saveasdir = os.path.join(config.plotsbasedir, "limits", foldername)
            try:
                os.makedirs(saveasdir)
            except OSError:
                pass
            plotlimits("higgsCombine_8TeV.MultiDimFit.mH125.root", "CMS_zz4l_fg4", os.path.join(saveasdir, "limit"))

if __name__ == "__main__":
    try:
        foldername = sys.argv[1]
    except IndexError:
        foldername = "test"
    runcombine(foldername)
