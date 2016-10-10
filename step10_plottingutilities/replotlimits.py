from helperstuff import config
from helperstuff.filemanager import cd, KeepWhileOpenFile
import os
import subprocess

def replotlimits(fileorfolder):
    if "bkp" in fileorfolder: return

    if os.path.isfile(fileorfolder) and fileorfolder.endswith(".txt"):
        with open(fileorfolder) as f:
            contents = f.readline().split()
            try:
                if contents[0] == "python" and contents[1] == "step9_runcombine.py":
                    print fileorfolder
                    filename = os.path.dirname(fileorfolder).replace("/", "")
                    with KeepWhileOpenFile(filename) as f:
                        print f, bool(f)
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
