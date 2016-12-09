#!/usr/bin/env python
from collections import OrderedDict
from copy import deepcopy
from helperstuff import config
from helperstuff.optimizesmoothing import TemplateIterate
from helperstuff.submitjob import submitjob
from helperstuff.templates import smoothingparametersfile, Template, templatesfiles
from helperstuff.utilities import cd, tfiles
import json
import os
import shutil
import subprocess
import sys
import ROOT

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
                shutil.move(templatesfile.templatesfile(), os.path.join(bkptemplatesdir, os.path.basename(templatesfile.templatesfile())))
            else:
                shutil.copy(templatesfile.templatesfile(), os.path.join(bkptemplatesdir, os.path.basename(templatesfile.templatesfile())))
        else:
            nchangedfiles += 1

        if os.path.exists(templatesfile.jsonfile()):
            shutil.move(templatesfile.jsonfile(), os.path.join(bkpjsondir, templatesfile.jsonfile()))

    jsonstring = json.dumps(smoothingparametersdict, sort_keys=True, indent=4, separators=(',', ': '))

    with open(smoothingparametersfile, "w") as f:
        f.write(jsonstring)

    Template.getsmoothingparametersdict(trycache=False)  #clear the cache
    tfiles.clear()

    print len(messages), nchangedfiles
    assert bool(messages) == bool(nchangedfiles)

    if messages:
        subprocess.check_call(["git", "add", smoothingparametersfile])
        message = ("Iterate smoothing parameters:\n\n"
                  + "\n".join("{}: {}".format(k, v) for k, v in messages.iteritems()))
        #print message
        #raw_input("press enter: ")
        subprocess.check_call(["git", "commit", "-m", message])
    else:
        subprocess.check_call(["git", "checkout", "--", smoothingparametersfile])

    return nchangedfiles

def runiteration(iternumber):
    iternumber = int(iternumber)

    smoothingparametersdict = Template.getsmoothingparametersdict(trycache=False)
    tfiles.clear()

    if "iteration" not in smoothingparametersdict:
        if iternumber != 0:
            raise IOError("Wrong iteration number {}!  smoothingparametersdict does not indicate an iteration number, so your iteration number should be 0.".format(iternumber))
        smoothingparametersdict["iteration"] = -1

    if iternumber == 0:
        if all(os.path.exists(_.templatesfile()) for _ in templatesfiles):
            return
        if any(os.path.exists(_.templatesfile()) or os.path.exists(_.jsonfile()) for _ in templatesfiles):
            raise IOError("For the 0th iteration, please delete all json and template files")

    if smoothingparametersdict["iteration"] == iternumber:
        #already iterated the values, just need to make the templates for this iteration
        nchangedfiles = 0
        for _ in templatesfiles:

            isgood = True
            if os.path.exists(_.templatesfile()):
                f = ROOT.TFile(_.templatesfile())
                if len(f.GetListOfKeys()) == 0:
                    isgood = False
                f.Close()
                if not isgood:
                    os.remove(_.templatesfile())

            if os.path.exists(_.templatesfile()+".tmp"):
                os.remove(_.templatesfile()+".tmp")

            if not os.path.exists(_.templatesfile()):
                nchangedfiles += 1

        if nchangedfiles == 0:
            return
    elif smoothingparametersdict["iteration"] == iternumber-1:
        if any(
               any(_2.smoothingparameters != [None, None, None] for _2 in _.templates())
                 and not os.path.exists(_.templatesfile()) for _ in templatesfiles
              ):
            raise IOError("Some templates have smoothing parameters defined but do not exist.  Please run the {}th iteration again.".format(iternumber-1))

        smoothingparametersdict["iteration"] = iternumber

        nchangedfiles = iterate(smoothingparametersdict, iternumber)
        converged = (nchangedfiles == 0)

        if converged:
            print
            print "***********************************" + "*"*(len(str(iternumber))-2)
            print "* converged after {} iterations!! *".format(iternumber)
            print "***********************************" + "*"*(len(str(iternumber))-2)
            print
            return True

    else:
        raise IOError("Wrong iteration number {}!  smoothingparametersdict says {}, so you should use either {} or {}".format(iternumber, smoothingparametersdict["iteration"], smoothingparametersdict["iteration"], smoothingparametersdict["iteration"]+1))

    with cd(config.repositorydir):
        subprocess.check_call(["python", "step4_makejson.py"])
        from helperstuff import ZX, CJLSTscripts #make sure the C++ stuff is compiled
                                                 #otherwise the jobs will try to compile them simultaneously
                                                 #and fail

        logdir = os.path.join(config.repositorydir, "step7_templates", "logs_iter{}".format(iternumber))
        try:
            os.makedirs(logdir)
        except OSError:
            pass
        outputfile = os.path.join(logdir, "log_{jobid}.out")

        ntries=3

        for ntry in range(ntries):
            print "submitting jobs..."
            jobids = []
            for i in range(nchangedfiles):
                jobids.append(submitjob("python step6_maketemplates.py", jobname="iter{}_job{}".format(iternumber, i), jobtime="10:0:0", outputfile=outputfile))

            submitjob("echo done", jobname="afteriteration", jobtime="0:0:10", interactive=True, waitids=jobids)

            print
            print "checking that the templates were all created..."

            allgood = True
            for templatesfile in templatesfiles:
                if not os.path.exists(templatesfile.templatesfile()):
                    allgood = False
                    print templatesfile.templatesfile(), "was not created."
                    if os.path.exists(templatesfile.templatesfile()+".tmp"):
                        os.remove(templatesfile.templatesfile()+".tmp")
                    continue

                f = ROOT.TFile(templatesfile.templatesfile())
                if len(f.GetListOfKeys()) == 0:
                    allgood = False
                    print templatesfile.templatesfile(), "contains nothing.  Removing it."
                    os.remove(templatesfile.templatesfile())
                f.Close()

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
        firstiter = previousruniter
    else:
        firstiter = 0

    for iternumber in range(firstiter, lastiter+1):
        if runiteration(iternumber):
            break
