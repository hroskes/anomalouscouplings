#!/usr/bin/env python

import json
import os

from helperstuff.templates import TemplatesFile, templatesfiles
from helperstuff.utilities import KeepWhileOpenFile, LSB_JOBID

def makejson(*args):
  templatesfile = TemplatesFile(*args)
  if templatesfile.copyfromothertemplatesfile is not None: return
  with KeepWhileOpenFile(templatesfile.jsonfile()+".tmp", message=LSB_JOBID()) as kwof:
    print templatesfile
    if not kwof: return
    if os.path.exists(templatesfile.jsonfile()): return

    jsonstring = json.dumps(templatesfile.getjson(), sort_keys=True, indent=4, separators=(',', ': '))
    filename = templatesfile.jsonfile()
    with open(filename, "w") as f:
      f.write(jsonstring)

if __name__ == "__main__":
    for templatesfile in templatesfiles:
        makejson(templatesfile)
