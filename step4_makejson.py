#!/usr/bin/env python
import argparse

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--overwrite", help="redo the json files that already exist", action="store_false", dest="keep")
    p.add_argument("--submitjobs", help="submit jobs", metavar="FILES_PER_JOB", type=int)
    p.add_argument("--filter", type=eval, default=lambda template: True)
    args = p.parse_args()

import datetime
import json
import os

from helperstuff import config
from helperstuff.submitjob import submitjob
from helperstuff.templates import TemplatesFile, templatesfiles
from helperstuff.utilities import KeepWhileOpenFile, LSB_JOBID, mkdir_p

def makejson(*args):
    templatesfile = TemplatesFile(*args)
    if templatesfile.copyfromothertemplatesfile is not None: return
    filename = templatesfile.jsonfile()
    mkdir_p(os.path.dirname(filename))
    with KeepWhileOpenFile(filename+".tmp") as f:
        if not f: return
        print templatesfile, datetime.datetime.now()

        jsonstring = json.dumps(templatesfile.getjson(), sort_keys=True, indent=4, separators=(',', ': '))
        with open(filename, "w") as f:
            f.write(jsonstring)

def submitjobs(filesperjob):
    i = 0
    for templatesfile in templatesfiles:
        filename = templatesfile.jsonfile()
        if os.path.exists(filename): continue
        if templatesfile.copyfromothertemplatesfile is not None: continue
        mkdir_p(os.path.dirname(filename))
        kwof = KeepWhileOpenFile(filename+".tmp")
        if not kwof.wouldbevalid:
            jobid = kwof.runningjobid
            if jobid: yield jobid
        if not i%filesperjob: yield submitjob("unbuffer "+os.path.join(config.repositorydir, "step4_makejson.py"), jobname="json"+str(i/filesperjob), jobtime="10:0:0", docd=True, queue="skylake")
        i += 1

if __name__ == "__main__":
    if args.submitjobs:
        list(submitjobs(args.submitjobs))
    else:
        for templatesfile in templatesfiles:
            if not args.filter(templatesfile): continue
            if args.keep and os.path.exists(templatesfile.jsonfile()): continue
            makejson(templatesfile)
