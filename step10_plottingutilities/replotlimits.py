import glob
from helperstuff import config
from helperstuff.filemanager import cd
import os
import subprocess

if __name__ == "__main__":
    with cd(config.repositorydir):
        for txtfile in glob.iglob(os.path.join(config.plotsbasedir, "limits", "*", "*.txt")):
            with open(txtfile) as f:
                contents = f.read().split()
                if contents[0] == "python" and contents[1] == "step9_runcombine.py":
                    subprocess.check_call(contents)
