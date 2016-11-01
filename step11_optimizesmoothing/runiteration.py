#!/usr/bin/env python
from collections import OrderedDict
from copy import deepcopy
from helperstuff import config
from helperstuff.filemanager import cd
from helperstuff.optimizesmoothing import TemplateIterate
from helperstuff.submitjob import submitjob
from helperstuff.templates import Template, templatesfiles
import json
import os
import shutil
import subprocess
import sys

def iterate(smoothingparametersdict, iternumber):
    bkptemplatesdir = os.path.join(config.repositorydir, "step7_templates", "bkp_iter{}".format(iternumber-1))
    bkpjsondir = os.path.join(config.repositorydir, "step5_json", "bkp_iter{}".format(iternumber-1))
    try:
        os.makedirs(bkptemplatesdir)
    except OSError:
        pass
    try:
        os.makedirs(bkpjsondir)
    except OSError:
        pass

    nchangedfiles = 0
    messages = OrderedDict()
    for templatesfile in templatesfiles:
        anychange = False
        for template in templatesfile.templates():

            template = TemplateIterate(template)
            nextiteration = template.getnextiteration()
            if nextiteration is None:
                continue
            nextiteration, message = nextiteration
            previousparameters = deepcopy(template.smoothingparameters)
            if nextiteration == previousparameters and os.path.exists(template.templatefile()):  #shouldn't be able to happen, but just in case
                continue
            anychange = True

            template.smoothingparameters[:] = nextiteration

            messages[template] = message

        if os.path.exists(templatesfile.templatesfile()):
            if anychange:
                nchangedfiles += 1
                shutil.move(templatesfile.templatesfile(), bkptemplatesdir)
            else:
                shutil.copy(templatesfile.templatesfile(), bkptemplatesdir)
        else:
            nchangedfiles += 1

        if os.path.exists(templatesfile.jsonfile()):
            shutil.move(templatesfile.jsonfile(), bkpjsondir)

    jsonstring = json.dumps(smoothingparametersdict, sort_keys=True, indent=4, separators=(',', ': '))

    with open(Template.smoothingparametersfile, "w") as f:
        f.write(jsonstring)

    Template.getsmoothingparametersdict(trycache=False)  #clear the cache

    print len(messages), nchangedfiles
    assert bool(messages) == bool(nchangedfiles)

    if messages:
        subprocess.check_call(["git", "add", Template.smoothingparametersfile])
        message = ("Iterate smoothing parameters:\n\n"
                  + "\n".join("{}: {}".format(k, v) for k, v in messages.iteritems()))
        subprocess.check_call(["git", "commit", "-m", message])

    return nchangedfiles

def runiteration(iternumber):
    iternumber = int(iternumber)

    smoothingparametersdict = Template.getsmoothingparametersdict(trycache=False)

    if "iteration" not in smoothingparametersdict:
        if iternumber != 0:
            raise IOError("Wrong iteration number {}!  smoothingparametersdict does not indicate an iteration number, so your iteration number should be 0.".format(iternumber))
        smoothingparametersdict["iteration"] = -1

    if smoothingparametersdict["iteration"] == iternumber:
        nchangedfiles = len(templatesfiles)
        pass #already iterated the values, just need to make the templates for this iteration
    elif smoothingparametersdict["iteration"] == iternumber-1:
        smoothingparametersdict["iteration"] = iternumber
 
        nchangedfiles = iterate(smoothingparametersdict, iternumber)
        converged = (nchangedfiles == 0)

        if converged:
            print
            print "***********************************" + "*"*len(str(iternumber))-2
            print "* converged after {} iterations!! *".format(iternumber)
            print "***********************************" + "*"*len(str(iternumber))-2
            print
            return True
    else:
        raise IOError("Wrong iteration number {}!  smoothingparametersdict says {}, so you should use either {} or {}".format(iternumber, smoothingparametersdict["iteration"], smoothingparametersdict["iteration"]+1))

    with cd(config.repositorydir):
        subprocess.check_call(["python", "step4_makejson.py"])
        from helperstuff import ZX, CJLSTscripts #make sure the C++ stuff is compiled
                                                 #otherwise the jobs will try to compile them simultaneously
                                                 #and fail

        ntries=3
        for ntry in range(ntries):
            print "submitting jobs..."
            jobids = []
            for i in range(nchangedfiles):
                jobids.append(submitjob("python step6_maketemplates.py", jobname="iter{}_job{}".format(iternumber, i), jobtime="2:0:0"))

            submitjob("echo done", jobname="afteriteration", jobtime="0:1:0", interactive=True, waitids=jobids)

            print
            print "checking that the templates were all created..."

            allgood = True
            for templatesfile in templatesfiles:
                if not os.path.exists(templatesfile.templatesfile()):
                    allgood = False
                    print templatesfile.templatesfile(), "was not created."
                    continue

                f = ROOT.TFile(templatesfile.templatesfile())
                if len(f.GetListOfKeys()) == 0:
                    allgood = False
                    print templatesfile.templatesfile(), "contains nothing.  Removing it."
                    os.remove(templatesfile.templatesfile())

            if allgood:
                print "...yes they were"
                break

            if ntry < ntries-1:
                print
                print "trying again"
                print

        else:
            raise RuntimeError("Tried to create templates", ntries, "times and failed.  Giving up.")

    return False

if __name__ == "__main__":
    lastiter = int(sys.argv[1])
    if "iteration" in Template.getsmoothingparametersdict():
        previousruniter = Template.getsmoothingparametersdict()["iteration"]
        if len(sys.argv) > 2:
            if sys.argv[2] == "remaketemplates":
                firstiter = previousruniter
            else:
                raise ValueError("unknown sys.argv[2] {}".format(sys.argv[2]))
        else:
           firstiter = previousruniter+1
    else:
        firstiter = 0

    for iternumber in range(firstiter, lastiter+1):
        if runiteration(iternumber):
            break
