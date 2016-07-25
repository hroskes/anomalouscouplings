import collections
from . import config
from .enums import EnumItem, Channel, DataTree, MultiEnum, MyEnum, Production, ProductionMode, Template
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
    enums = [ProductionMode, Channel, Luminosity, Production]
    def getrate(self):
        if self.productionmode == "ZX" and self.production == "160720":
            if self.channel == "4e":    return 1.36 * float(self.luminosity)/7.65
            if self.channel == "4mu":   return 1.64 * float(self.luminosity)/7.65
            if self.channel == "2e2mu": return 2.81 * float(self.luminosity)/7.65
        if self.productionmode == "ZX" and self.production == "160225": #email from Simon, Feb 9 at 4:56 PM, "inputs for the cards"
            if self.channel == "4e":    return (0.311745 + 0.0106453) * float(self.luminosity)/2.8
            if self.channel == "4mu":   return 0.408547 * float(self.luminosity)/2.8
            if self.channel == "2e2mu": return (0.716686 + 0.0199815) * float(self.luminosity)/2.8

        if self.productionmode == "ggH":
            result = Template("fa3", self.productionmode, self.channel, "0+", self.production).gettemplate().Integral()*float(self.luminosity)
            for productionmode in "VBF", "WplusH", "WminusH", "ZH", "ttH":
                sample = Sample(productionmode, "0+", self.production)
                f = tfiles[sample.withdiscriminantsfile()]
                t = f.candTree
                ZZFlav = self.channel.ZZFlav
                additionalxsec = 0
                for event in t:
                    if config.m4lmin < t.ZZMass < config.m4lmax and t.Z1Flav*t.Z2Flav == ZZFlav:
                        additionalxsec += getattr(t, sample.weightname())
                result += additionalxsec * float(self.luminosity)
            return result

        result = Template("fa3", self.productionmode, self.channel, self.production).gettemplate().Integral()*float(self.luminosity)
        return result

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
    return Template(*args).gettemplate()

def getdatatree(*args):
    return tfiles[DataTree(*args).treefile].candTree
    #it's empty if we don't actually unblind

def discriminantnames(*args):
    theset = set()
    for production in config.productionsforcombine:
        result = Template("ggH", "0+", "2e2mu", production, *args).discriminants()
        theset.add(result)
    assert len(theset) == 1  #sanity check
    return result
