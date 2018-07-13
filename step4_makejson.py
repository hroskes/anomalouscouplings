#!/usr/bin/env python
import argparse

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--overwrite", help="keep the json files that already exist", action="store_false", dest="keep")
    p.add_argument("--submitjobs", help="submit jobs", metavar="FILES_PER_JOB", type=int)
    args = p.parse_args()

import datetime
import json
import os

from helperstuff import config
from helperstuff.submitjob import submitjob
from helperstuff.templates import TemplatesFile, templatesfiles
from helperstuff.utilities import KeepWhileOpenFile, LSB_JOBID

def makejson(*args):
    templatesfile = TemplatesFile(*args)
    if templatesfile.copyfromothertemplatesfile is not None: return
    with KeepWhileOpenFile(templatesfile.jsonfile()+".tmp") as f:
        if not f: return
        print templatesfile, datetime.datetime.now()

        jsonstring = json.dumps(templatesfile.getjson(), sort_keys=True, indent=4, separators=(',', ': '))
        filename = templatesfile.jsonfile()
        with open(filename, "w") as f:
            f.write(jsonstring)

def submitjobs(filesperjob):
    i = 0
    for templatesfile in templatesfiles:
        if os.path.exists(templatesfile.jsonfile()) or not KeepWhileOpenFile(templatesfile.jsonfile()+".tmp").wouldbevalid: continue
        if not i%filesperjob: yield submitjob("unbuffer "+os.path.join(config.repositorydir, "step4_makejson.py"), jobname="json"+str(i/filesperjob), jobtime="10:0:0", docd=True)
        i += 1

if __name__ == "__main__":
    if args.submitjobs:
        list(submitjobs(args.submitjobs))
    else:
        for templatesfile in templatesfiles:
            if args.keep and os.path.exists(templatesfile.jsonfile()): continue
            makejson(templatesfile)
