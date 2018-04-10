#!/usr/bin/env python
from helperstuff import config
from helperstuff.utilities import cd, KeepWhileOpenFile
import os
import subprocess

def replotlimits(fileorfolder):
    if "bkp" in fileorfolder.replace(config.plotsbasedir, ""): return

    print fileorfolder
    if os.path.isfile(fileorfolder) and fileorfolder.endswith(".txt"):
        if "February" not in fileorfolder: return
        with open(fileorfolder) as f:
            contents = f.readline().split()
            try:
                if contents[0] == "python" and contents[1] in ["step9_runcombine.py", "./step9_runcombine.py"]:
                    print fileorfolder
                    filename = os.path.dirname(fileorfolder).replace("/", "")
                    with KeepWhileOpenFile(filename) as f:
                        if f:
                            subprocess.check_call(contents)
                if contents[0] == "python" and contents[1] in ["mergeplots.py", "./mergeplots.py"]:
                    print fileorfolder
                    filename = os.path.dirname(fileorfolder).replace("/", "")
                    with KeepWhileOpenFile(filename) as f, cd("step10_plottingutilities"):
                        if f:
                            subprocess.check_call(contents)
            except IndexError:
                print fileorfolder
                raise

    if os.path.isdir(fileorfolder):
        for filename in os.listdir(fileorfolder):
            replotlimits(os.path.join(fileorfolder, filename))

if __name__ == "__main__":
    with cd(config.repositorydir):
        replotlimits(os.path.join(config.plotsbasedir, "limits"))
