#!/usr/bin/env python
from helperstuff import utilities
import os
import shutil
import subprocess
import tempfile

print """
Please set this repository's location and the directory to store plots in helperstuff/config.py
Checking to see if you've done it already..."""
import helperstuff.config
if os.path.realpath(os.getcwd()) != os.path.realpath(helperstuff.config.repositorydir):
    raise OSError("You have a setup at {}, but you're not there!".format(os.path.realpath(helperstuff.config.repositorydir)))
print """Yes, config is set up!
"""

print "Initiating git submodules..."
subprocess.check_call(["git", "submodule", "update", "--init", "--recursive"])

if set(os.listdir("CMSSW_10_2_5")) != set(["src", ".gitignore"]):
    print """CMSSW_10_2_5 area already set up."""
    print
else:
    print """Now setting up CMSSW_10_2_5 area for combine..."""
    #move files out of CMSSW_10_2_5
    tmpdir = "tmpCMSSW_10_2_5"
    utilities.mkdir_p(tmpdir)
    shutil.move("CMSSW_10_2_5/src", os.path.join(tmpdir, "src"))
    shutil.move("CMSSW_10_2_5/.gitignore", os.path.join(tmpdir, ".gitignore"))
    try:
        os.rmdir("CMSSW_10_2_5")
        os.environ["SCRAM_ARCH"] = "slc{}_amd64_gcc630".format(helperstuff.config.slcversion)
        subprocess.check_call(["scram", "p", "CMSSW", "CMSSW_10_2_5"])
    finally:
        if os.path.exists("CMSSW_10_2_5/src"):
            os.rmdir("CMSSW_10_2_5/src")
        shutil.move(os.path.join(tmpdir, "src"), "CMSSW_10_2_5/src")
        shutil.move(os.path.join(tmpdir, ".gitignore"), "CMSSW_10_2_5/.gitignore")
        os.rmdir(tmpdir)

    print """CMSSW area is set up"""

print
print "Compiling MELA..."

with utilities.cd("CMSSW_10_2_5/src/ZZMatrixElement"):
    with open(".git") as f:
        gitdir = f.read().replace("gitdir: ", "").strip()
    with open(os.path.join(gitdir, "info", "sparse_checkout"), "w") as f:
        f.write("/MELA\n/setup.sh\n/README.md\n/.gitignore")
    subprocess.check_call(["git", "config", "core.sparsecheckout", "true"])
    subprocess.check_call(["git", "read-tree", "-mu", "HEAD"])
    subprocess.check_call("eval $(scram ru -sh) && ./setup.sh -j 10", shell=True)

print
print """Compiling CMSSW..."""

with utilities.cd("CMSSW_10_2_5/src"):
    subprocess.check_call(["scram", "b", "-j", "10"])

print "Compiling NIS_summary..."
subprocess.check_call("cd CMSSW_10_2_5 && eval $(scram ru -sh) && cd ../step10_plottingutilities/NIS_summary && make", shell=True)

print """Checking that python dependencies are installed..."""
try:
  import uncertainties
except ImportError:
  print "Installing uncertainties..."
  with utilities.cd("CMSSW_10_2_5"):
    subprocess.check_call("eval $(scram ru -sh) && pip install --user uncertainties", shell=True)
