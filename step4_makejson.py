from Alignment.OfflineValidation.TkAlAllInOneTool.helperFunctions import replaceByMap  #easiest place to get it
from helperstuff.samples import Sample
from helperstuff.enums import Channel, channels, treesystematics, TemplatesFile
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
    if templatesfile.signalorbkg == "sig":
        samplesdict = {
                       Sample("ggH", "0+"): {
                                             "templatename": "0Plus",
                                             "weightname": "MC_weight_ggH_g1",
                                             "scalefactor": str(1./6),
                                            },
                       Sample("ggH", "0-"): {
                                             "templatename": "0Minus",
                                             "weightname": "MC_weight_ggH_g4",
                                             "scalefactor": str(1./6),
                                            },
                       Sample("ggH",
                                 "fa30.5"): {
                                             "templatename": "g1g4",
                                             "weightname": "MC_weight_ggH_g1g4",
                                             "scalefactor": str(2./6),
                                            },
                  }
    else:
        samplesdict = {
                       Sample("qqZZ"):      {
                                             "templatename": "qqZZ",
                                             "weightname": "MC_weight_qqZZ",
                                             "scalefactor": "1",
                                            },
                       Sample("ggZZ",
                                  "2e2mu"): {
                                             "templatename": "ggZZ",
                                             "weightname": "MC_weight_ggZZ",
                                             "scalefactor": "1",
                                            },
                       Sample("ZX"):        {
                                             "templatename": "ZX",
                                             "weightname": "MC_weight_ZX",
                                             "scalefactor": "1",
                                            },
                  }
    domirror = {
                Sample("ggH", "0+"): True,
                Sample("ggH", "0-"): True,
                Sample("ggH", "fa30.5"): False,
                Sample("ggZZ", "2e2mu"): True,
                Sample("qqZZ"): True,
                Sample("ZX"): True,
               }
    flavorproducts = {
                      Channel("2e2mu"): str(13*13*11*11),
                      Channel("4e"): str(11*11*11*11),
                      Channel("4mu"): str(13*13*13*13),
                     }
    generalmap = {
                  "flavor": str(templatesfile.channel),
                  "flavorproduct": flavorproducts[templatesfile.channel],
                  "systematic": templatesfile.systematic.appendname(),
                 }


    for sample, samplemap in samplesdict.iteritems():
        if templatesfile.signalorbkg == "sig":
            sampleslist = [
                           Sample("ggH", "0+"),
                           Sample("ggH", "a2"),
                           Sample("ggH", "0-"),
                           Sample("ggH", "L1"),
                           Sample("ggH", "fa20.5"),
                           Sample("ggH", "fa30.5"),
                           #Sample("ggH", "fL10.5"),   #NOT fL1 for now
                          ]
        elif sample.productionmode == "ggZZ":
            sampleslist = [
                           Sample("ggZZ", "4e"),
                           Sample("ggZZ", "2e2mu"),
                           Sample("ggZZ", "2e2tau"),
                           Sample("ggZZ", "4mu"),
                           Sample("ggZZ", "2mu2tau"),
                           Sample("ggZZ", "4tau"),
                          ]
        else:
            sampleslist = [sample]

        fileslist = [os.path.basename(s.withdiscriminantsfile()) for s in sampleslist]

        repmap = samplemap
        repmap.update(generalmap)

        with open("jsontemplates/basetemplate.json") as f:
            basetemplate = f.read()
        basetemplate = replaceByMap(basetemplate, repmap)
        basejson = jsonloads(basetemplate)
        basejson["templates"][0]["files"] = fileslist

        if domirror[sample]:
            with open("jsontemplates/mirror.json") as f:
                mirrortemplate = f.read()
            mirrortemplate = replaceByMap(mirrortemplate, repmap)
            mirrorjson = jsonloads(mirrortemplate)
            basejson["templates"] += mirrorjson["templates"]
            #flooring is done to the mirrored template
        else:
            basejson["templates"][0]["postprocessing"].append({"type": "floor"})

        if templatesfile.signalorbkg == "bkg":
            if sample == Sample("ZX"):
                basejson["templates"][0]["postprocessing"][0]["entriesperbin"] = 5
            else:
                basejson["templates"][0]["postprocessing"][0]["entriesperbin"] = 20
            #del basejson["templates"][0]["postprocessing"][0]

        jsondict["templates"] += basejson["templates"]

    if templatesfile.signalorbkg == "sig":
        with open("jsontemplates/int.json") as f:
            inttemplate = f.read()
        intjson = jsonloads(inttemplate)
        jsondict["templates"] += intjson["templates"]

    jsonstring = json.dumps(jsondict, sort_keys=True, indent=4, separators=(',', ': '))
    jsonstring = replaceByMap(jsonstring, generalmap)
    filename = templatesfile.jsonfile()
    with open(filename, "w") as f:
        f.write(jsonstring)

if __name__ == "__main__":
    for channel in channels:
        for systematic in treesystematics:
            makejson(channel, systematic, "signal")
        makejson(channel, "bkg")
