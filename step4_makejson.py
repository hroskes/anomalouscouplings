#!/usr/bin/env python
from helperstuff.templates import TemplatesFile, templatesfiles
import json

def makejson(*args):
    templatesfile = TemplatesFile(*args)
    if templatesfile.copyfromothertemplatesfile is not None: return
    print templatesfile

    jsonstring = json.dumps(templatesfile.getjson(), sort_keys=True, indent=4, separators=(',', ': '))
    filename = templatesfile.jsonfile()
    with open(filename, "w") as f:
        f.write(jsonstring)

if __name__ == "__main__":
    for templatesfile in templatesfiles:
        makejson(templatesfile)
