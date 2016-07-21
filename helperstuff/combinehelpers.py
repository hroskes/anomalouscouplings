import collections
from . import config
from .enums import EnumItem, Channel, MultiEnum, MyEnum, Production, ProductionMode, Template
from .filemanager import tfiles
from .samples import Sample
import ROOT

class LuminosityType(MyEnum):
    enumname = "luminositytype"
    enumitems = (
                 EnumItem("fordata"),
                 EnumItem("forexpectedscan"),
                )

class Luminosity(MultiEnum):
    enumname = "luminosity"
    enums = [LuminosityType, Production]
    def check(self, *args):
        super(Luminosity, self).check(self, *args, dontcheck=[Production])
        if self.luminositytype == "fordata":
            super(Luminosity, self).check(self, *args)
    def __float__(self):
        if self.luminositytype == "forexpectedscan": return float(config.expectedscanluminosity)
        if self.luminositytype == "fordata": return float(self.production.dataluminosity)

class __Rate(MultiEnum):
    enums = [ProductionMode, Channel, Luminosity]
    def getrate(self):
        if self.productionmode == "ggZZ":
            if self.channel == "4e":    return 0.40 * float(self.luminosity)/10
            if self.channel == "4mu":   return 0.77 * float(self.luminosity)/10
            if self.channel == "2e2mu": return 0.66 * float(self.luminosity)/10
        if self.productionmode == "qqZZ":
            if self.channel == "4e":    return 3.26 * float(self.luminosity)/10
            if self.channel == "4mu":   return 7.17 * float(self.luminosity)/10
            if self.channel == "2e2mu": return 8.77 * float(self.luminosity)/10
        if self.productionmode == "ZX":
            if self.channel == "4e":    return 1.36 * float(self.luminosity)/7.6
            if self.channel == "4mu":   return 1.64 * float(self.luminosity)/7.6
            if self.channel == "2e2mu": return 2.81 * float(self.luminosity)/7.6

        if self.productionmode == "ggH":
            result = Template("fa3", self.productionmode, self.channel, "0+", config.productionforcombine).gettemplate().Integral()*float(self.luminosity)
            for productionmode in "VBF", "WplusH", "WminusH", "ZH", "ttH":
                sample = Sample(productionmode, "0+", config.productionforcombine)
                f = tfiles[sample.withdiscriminantsfile()]
                t = f.candTree
                ZZFlav = self.channel.ZZFlav
                additionalxsec = 0
                for event in t:
                    if 105 < t.ZZMass < 140 and t.Z1Flav*t.Z2Flav == ZZFlav:
                        additionalxsec += getattr(t, sample.weightname())
                result += additionalxsec * float(self.luminosity)
            return result

        assert False
        return Template("fa3", self.productionmode, self.channel, config.productionforcombine).gettemplate().Integral()*float(self.luminosity)

    def __float__(self):
        return self.getrate()

__cache = {}
def getrate(*args):
    rate = __Rate(*args)
    if rate not in __cache:
        __cache[rate] = float(rate)
    return __cache[rate]

def getrates(*args):
    ggH, qqZZ, ggZZ, ZX = getrate("ggH", *args), getrate("qqZZ", *args), getrate("ggZZ", *args), getrate("ZX", *args)

    result =  "rate {} {} {} {}".format(ggH, qqZZ, ggZZ, ZX)
    return result

def gettemplate(*args):
    return Template(config.productionforcombine, *args).gettemplate()

def getdatatree(channel):
    channel = Channel(channel)
    return tfiles[Sample("data", config.productionforcombine, "unblind").withdiscriminantsfile()].candTree
    #unblind is empty if we don't actually unblind
