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

if set(os.listdir("CMSSW_8_1_0")) != set(["src", ".gitignore"]):
    print """CMSSW_8_1_0 area already set up."""
    print
else:
    print """Now setting up CMSSW_8_1_0 area for combine..."""
    #move files out of CMSSW_8_1_0
    tmpdir = "tmpCMSSW_8_1_0"
    utilities.mkdir_p(tmpdir)
    shutil.move("CMSSW_8_1_0/src", os.path.join(tmpdir, "src"))
    shutil.move("CMSSW_8_1_0/.gitignore", os.path.join(tmpdir, ".gitignore"))
    try:
        os.rmdir("CMSSW_8_1_0")
        os.environ["SCRAM_ARCH"] = "slc6_amd64_gcc493"
        subprocess.check_call(["scram", "p", "CMSSW", "CMSSW_8_1_0"])
    finally:
        if os.path.exists("CMSSW_8_1_0/src"):
            os.rmdir("CMSSW_8_1_0/src")
        shutil.move(os.path.join(tmpdir, "src"), "CMSSW_8_1_0/src")
        shutil.move(os.path.join(tmpdir, ".gitignore"), "CMSSW_8_1_0/.gitignore")

    print """CMSSW area is set up"""
    print

print """Compiling CMSSW..."""

with utilities.cd("CMSSW_8_1_0/src"):
    subprocess.check_call(["scram", "b", "-j", "10"])

print "Compiling TemplateBuilder..."

subprocess.check_call("cd CMSSW_8_1_0 && eval $(scram ru -sh) && cd ../TemplateBuilder && make", shell=True)
gitignore = """
obj/*
buildTemplate.exe
.gitignore
"""
with open("TemplateBuilder/.gitignore", "w") as f:
    f.write(gitignore)

print "Compiling NIS_summary..."
subprocess.check_call("cd CMSSW_8_1_0 && eval $(scram ru -sh) && cd ../step10_plottingutilities/NIS_summary && make", shell=True)
