import collections
from . import config
from .enums import Channel, MultiEnum, ProductionMode, Template
from .filemanager import tfiles
from .samples import Sample
import ROOT


class __Rate(MultiEnum):
    enums = (ProductionMode, Channel)
    def getrate(self):
        if self.productionmode == "ggZZ":
            if self.channel == "4e":    return 0.40 * config.luminosity/10
            if self.channel == "4mu":   return 0.77 * config.luminosity/10
            if self.channel == "2e2mu": return 0.66 * config.luminosity/10
        if self.productionmode == "qqZZ":
            if self.channel == "4e":    return 3.26 * config.luminosity/10
            if self.channel == "4mu":   return 7.17 * config.luminosity/10
            if self.channel == "2e2mu": return 8.77 * config.luminosity/10
        if self.productionmode == "ZX":
            if self.channel == "4e":    return 2.196 * config.luminosity/10
            if self.channel == "4mu":   return 3.003 * config.luminosity/10
            if self.channel == "2e2mu": return 3.116 * config.luminosity/10

        if self.productionmode == "ggH":
            result = Template("fa3", self, "0+", config.productionforsignalrates).gettemplate().Integral()*config.luminosity
            for productionmode in "VBF", "WplusH", "WminusH", "ZH", "ttH":
                sample = Sample(productionmode, "0+", config.productionforsignalrates)
                f = tfiles[sample.withdiscriminantsfile()]
                t = f.candTree
                ZZFlav = self.channel.ZZFlav()
                additionalxsec = 0
                for event in t:
                    if 105 < t.ZZMass < 140 and t.Z1Flav*t.Z2Flav == ZZFlav:
                        additionalxsec += getattr(t, sample.weightname())
                result += additionalxsec * config.luminosity
            return result

        assert False
        return Template("fa3", self, config.productionforsignalrates).gettemplate().Integral()*config.luminosity

    def __float__(self):
        return self.getrate()

__cache = {}
def getrate(*args):
    rate = __Rate(*args)
    if rate not in __cache:
        __cache[rate] = float(rate)
    return __cache[rate]

def getrates(flavor):
    ggH, qqZZ, ggZZ, ZX = getrate("ggH", flavor), getrate("qqZZ", flavor), getrate("ggZZ", flavor), getrate("ZX", flavor)

    result =  "rate {} {} {} {}".format(ggH, qqZZ, ggZZ, ZX)
    return result

def gettemplate(*args):
    return Template(config.productionforcombine, *args).gettemplate()

def getdatatree(channel):
    channel = Channel(channel)
    return tfiles[Sample("data", config.productionforcombine, "unblind").withdiscriminantsfile()].candTree
    #unblind is empty if we don't actually unblind
