#!/usr/bin/env python
import argparse

p = argparse.ArgumentParser()
p.add_argument("--keep", help="keep the json files that already exist", action="store_true")
args = p.parse_args()

import datetime
import json
import os

from helperstuff.templates import TemplatesFile, templatesfiles
from helperstuff.utilities import KeepWhileOpenFile

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

if __name__ == "__main__":
    for templatesfile in templatesfiles:
        if templatesfile.production == "180416_Ulascan": continue
        if args.keep and os.path.exists(templatesfile.jsonfile()): continue
        makejson(templatesfile)
