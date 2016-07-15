from helperstuff import constants
from helperstuff.samples import Sample
from helperstuff.enums import analyses, Channel, channels, productions, treesystematics, TemplatesFile
import json
import os

def jsonloads(jsonstring):
    try:
        return json.loads(jsonstring)
    except:
        print jsonstring
        raise

def makejson(*args):
    templatesfile = TemplatesFile(*args)
    jsondict = {
                "inputDirectory": "step3_withdiscriminants/",
                "outputFile": templatesfile.templatesfile(),
                "templates": [],
               }

    for template in templatesfile.templates():
        jsondict["templates"] += template.getjson()["templates"]

    jsonstring = json.dumps(jsondict, sort_keys=True, indent=4, separators=(',', ': '))
    filename = templatesfile.jsonfile()
    with open(filename, "w") as f:
        f.write(jsonstring)

if __name__ == "__main__":
    for production in productions:
        for analysis in analyses:
            for channel in channels:
                for systematic in treesystematics:
                    print channel, analysis, systematic
                    makejson(channel, systematic, "signal", analysis, production)
                makejson(channel, "bkg", analysis, production)
