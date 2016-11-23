#!/usr/bin/env python
from helperstuff import filemanager
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

if set(os.listdir("CMSSW_7_6_5")) != set(["src", ".gitignore"]):
    print """CMSSW_7_6_5 area already set up."""
    print
else:
    print """Now setting up CMSSW_7_6_5 area for combine..."""
    #move files out of CMSSW_7_6_5
    tmpdir = tempfile.mkdtemp()
    shutil.move("CMSSW_7_6_5/src", os.path.join(tmpdir, "src"))
    shutil.move("CMSSW_7_6_5/.gitignore", os.path.join(tmpdir, ".gitignore"))
    try:
        os.rmdir("CMSSW_7_6_5")
        os.environ["SCRAM_ARCH"] = "slc6_amd64_gcc493"
        subprocess.check_call(["scram", "p", "CMSSW", "CMSSW_7_6_5"])
    finally:
        if os.path.exists("CMSSW_7_6_5/src"):
            os.rmdir("CMSSW_7_6_5/src")
        shutil.move(os.path.join(tmpdir, "src"), "CMSSW_7_6_5/src")
        shutil.move(os.path.join(tmpdir, ".gitignore"), "CMSSW_7_6_5/.gitignore")

    with filemanager.cd("CMSSW_7_6_5/src"):
        subprocess.check_call(["scram", "b"])

    print """CMSSW area is set up"""

print "Compiling TemplateBuilder"

subprocess.check_call("cd CMSSW_7_6_5 && eval $(scram ru -sh) && cd ../TemplateBuilder && make", shell=True)
gitignore = """
obj/*
buildTemplate.exe
.gitignore
"""
with open("TemplateBuilder/.gitignore", "w") as f:
    f.write(gitignore)
