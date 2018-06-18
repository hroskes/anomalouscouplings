#!/usr/bin/env python
import argparse

p = argparse.ArgumentParser()
p.add_argument("--keep", help="keep the json files that already exist", action="store_true")
p.add_argument("--submitjobs", help="submit jobs", action="store_true")
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
    with KeepWhileOpenFile(templatesfile.jsonfile()+".tmp", message=LSB_JOBID()) as f:
        if not f: return
        print templatesfile, datetime.datetime.now()

        jsonstring = json.dumps(templatesfile.getjson(), sort_keys=True, indent=4, separators=(',', ': '))
        filename = templatesfile.jsonfile()
        with open(filename, "w") as f:
            f.write(jsonstring)

def submitjobs():
    i = 0
    for templatesfile in templatesfiles:
        if os.path.exists(templatesfile.jsonfile()) or os.path.exists(templatesfile.jsonfile()+".tmp"): continue
        submitjob("unbuffer "+os.path.join(config.repositorydir, "step4_makejson.py")+" --keep", jobname="json"+str(i), jobtime="30:0:0", docd=True)
        i += 1

if __name__ == "__main__":
    if args.submitjobs:
        submitjobs()
    else:
        for templatesfile in templatesfiles:
            if args.keep and os.path.exists(templatesfile.jsonfile()): continue
            makejson(templatesfile)
