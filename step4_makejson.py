from Alignment.OfflineValidation.TkAlAllInOneTool.helperFunctions import replaceByMap  #easiest place to get it
from helperstuff.samples import Sample
from helperstuff.enums import Channel
import json

def jsonloads(jsonstring):
    try:
        return json.loads(jsonstring)
    except:
        print jsonstring
        raise

samples = [
           Sample("ggH", "0+"),
           Sample("ggH", "a2"),
           Sample("ggH", "0-"),
           Sample("ggH", "L1"),
           Sample("ggH", "fa20.5"),
           Sample("ggH", "fa30.5"),
           #Sample("ggH", "fL10.5"),   #NOT L1 for now
          ]

def makejson(flavor):
    flavor = Channel(flavor)
    jsondict = {
                "inputDirectory": "step3_withdiscriminants/",
                "outputFile": "step5_templates/.oO[flavor]Oo._fa3Adap_new.root",
                "templates": []
               }
    samplesdict = {
                   Sample("ggH", "0+"): {
                                         "templatename": "0Plus",
                                         "weightname": "MC_weight_g1",
                                         "scalefactor": "1",
                                        },
                   Sample("ggH", "0-"): {
                                         "templatename": "0Minus",
                                         "weightname": "MC_weight_g4",
                                         "scalefactor": "1",
                                        },
                   Sample("ggH",
                             "fa30.5"): {
                                         "templatename": "g1g4",
                                         "weightname": "MC_weight_g1g4",
                                         "scalefactor": "2",
                                        },
                  }
    domirror = {
                Sample("ggH", "0+"): True,
                Sample("ggH", "0-"): True,
                Sample("ggH", "fa30.5"): False,
               }
    flavorproducts = {
                      Channel("2e2mu"): str(13*13*11*11),
                      Channel("4e"): str(11*11*11*11),
                      Channel("4mu"): str(13*13*13*13),
                     }
    flavormap = {
                 "flavor": str(flavor),
                 "flavorproduct": flavorproducts[flavor],
                }

    fileslist = [sample.withdiscriminantsfile() for sample in samples]

    for sample, samplemap in samplesdict.iteritems():
        repmap = samplemap
        repmap.update(flavormap)

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


        jsondict["templates"] += basejson["templates"]

    with open("jsontemplates/int.json") as f:
        inttemplate = f.read()
    intjson = jsonloads(inttemplate)
    jsondict["templates"] += intjson["templates"]

    jsonstring = json.dumps(jsondict, sort_keys=True, indent=4, separators=(',', ': '))
    jsonstring = replaceByMap(jsonstring, repmap)
    with open("step5_json/templates_{}.json".format(flavor), "w") as f:
        f.write(jsonstring)

if __name__ == "__main__":
    makejson("2e2mu")
    makejson("4e")
    makejson("4mu")
