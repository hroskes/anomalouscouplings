import collections
from helperstuff import config
from helperstuff.enums import Channel, channels, Template
from helperstuff.filemanager import tfiles
from helperstuff.samples import Sample
import ROOT

def getrates(flavor, analysis):
    flavor = Channel(flavor)
    ggH = Template(flavor, analysis, "ggH", "0+").gettemplate().Integral()*config.luminosity
    #other signal samples
    for productionmode in "VBF", "WplusH", "WminusH", "ZH", "ttH":
        sample = Sample(productionmode, "0+")
        f = tfiles[sample.withdiscriminantsfile()]
        t = f.candTree
        ZZFlav = flavor.ZZFlav()
        additionalxsec = 0
        for event in t:
            if 105 < t.ZZMass < 140 and t.Z1Flav*t.Z2Flav == ZZFlav:
                additionalxsec += getattr(t, sample.weightname())
        ggH += additionalxsec * config.luminosity


    qqZZ = Template(flavor, analysis, "qqZZ").gettemplate().Integral()*config.luminosity
    ggZZ = Template(flavor, analysis, "ggZZ").gettemplate().Integral()*config.luminosity
    ZX = Template(flavor, analysis, "ZX").gettemplate().Integral() * (
         #no idea about the absolute
         #rescale to Moriond
             (0.408547+0.311745+0.0106453+0.716686+0.0199815)  #email from Simon, "Inputs for the cards", Feb 9 at 4:56AM
              / sum(Template(c, "ZX", analysis).gettemplate().Integral() for c in channels)
         # * ratio of luminosities
             * config.luminosity / 2.8
         )

    ######################
    #       PATCH        #
    ######################
    #qqZZ and ggZZ are from 76X MC
    #7X from https://indico.cern.ch/event/515361/contributions/2223587/attachments/1302276/1944656/20160701_hzz4l_zxSS_v1.pdf
    if flavor == "4e":    qqZZ, ggZZ, ZX = 3.26, 0.40, 2.196
    if flavor == "4mu":   qqZZ, ggZZ, ZX = 7.17, 0.77, 3.003
    if flavor == "2e2mu": qqZZ, ggZZ, ZX = 8.77, 0.66, 3.116

    result =  "rate {} {} {} {}".format(ggH, qqZZ, ggZZ, ZX)
    return result

def gettemplate(*args):
    return Template(*args).gettemplate()

def getdatatree(channel):
    channel = Channel(channel)
    return tfiles[Sample("data").withdiscriminantsfile()].candTree
