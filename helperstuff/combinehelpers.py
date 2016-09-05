import collections
import config
from enums import Analysis, categories, Category, Channel, EnumItem, MultiEnum, MyEnum, Production, ProductionMode
from filemanager import tfiles
import os
import ROOT
from samples import ReweightingSample, Sample
from templates import DataTree, IntTemplate, Template
import yaml

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
    enums = [ProductionMode, Channel, Luminosity, Production, Category]

    def getrate(self):
        if self.productionmode in ["ggH", "VBF"]:
            hypothesis = "0+"
        else:
            hypothesis = None

        categorysum = 0
        for category in categories:
            categorysum += Template("fa3", "prod+dec", "D_int_decay" if category == "VBF2jTaggedIchep16" else None, hypothesis, category, self.productionmode, self.channel, self.production).gettemplate().Integral()*float(self.luminosity)

        myintegral = Template("fa3", "prod+dec", "D_int_decay" if self.category == "VBF2jTaggedIchep16" else None, hypothesis, self.category, self.productionmode, self.channel, self.production).gettemplate().Integral()*float(self.luminosity)

        return self.totalrate() * myintegral/categorysum

    def totalrate(self):
        if self.productionmode != "VBF bkg":
            return self.yamlrate()
        if self.productionmode == "ZX" and self.production in ("160725", "160729"):
            if self.channel == "4e":    return 2.39 * float(self.luminosity)/12.9
            if self.channel == "4mu":   return 3.66 * float(self.luminosity)/12.9
            if self.channel == "2e2mu": return 6.29 * float(self.luminosity)/12.9
        if self.productionmode == "ZX" and self.production == "160720":
            if self.channel == "4e":    return 1.36 * float(self.luminosity)/7.65
            if self.channel == "4mu":   return 1.64 * float(self.luminosity)/7.65
            if self.channel == "2e2mu": return 2.81 * float(self.luminosity)/7.65
        if self.productionmode == "ZX" and self.production == "160225": #email from Simon, Feb 9 at 4:56 PM, "inputs for the cards"
            if self.channel == "4e":    return (0.311745 + 0.0106453) * float(self.luminosity)/2.8
            if self.channel == "4mu":   return 0.408547 * float(self.luminosity)/2.8
            if self.channel == "2e2mu": return (0.716686 + 0.0199815) * float(self.luminosity)/2.8
        if self.productionmode == "ZX":
            assert False

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

        result = 0
        for category in categories:
            result += Template("fa3", "prod+dec", "D_int_decay" if category == "VBF2jTaggedIchep16" else None, category, self.productionmode, self.channel, self.production).gettemplate().Integral()*float(self.luminosity)
        return result

    def yamlrate(self):
        if self.production.year == 2016:
            filename = os.path.join(config.repositorydir, "helperstuff", "Datacards13TeV_ICHEP2016", "LegoCards", "configs", "inputs", "yields_per_tag_category_13TeV_{}.yaml".format(self.channel))
            tags = ["UnTagged", "VBF1jTagged", "VBF2jTagged", "VHLeptTagged", "VHHadrTagged", "ttHTagged"]
        elif self.production.year == 2015:
            filename = os.path.join(config.repositorydir, "helperstuff", "Datacards13TeV_Moriond2016", "LegoCards", "configs", "inputs", "yields_per_tag_category_13TeV_{}.yaml".format(self.channel))
            tags = ["UnTagged", "VBFTagged"]
        with open(filename) as f:
            y = yaml.load(f)
        with open(filename) as f:
            for line in f:
                if "fb-1" in line:
                    lumi = float(line.split("=")[1].split("fb-1")[0])
                    break
            else:
                raise IOError("No luminosity in {}".format(filename))
        if self.productionmode == "ggH":
            productionmodes = ["ggH", "qqH", "WH", "ZH", "ttH"]
        elif self.productionmode == "ZX":
            productionmodes = ["zjets"]
        else:
            productionmodes = [str(self.productionmode)]
        rate = 0

        for tag in tags:
            for p in productionmodes:
                try:
                    rate += float(y[tag][p]) * float(self.luminosity) / lumi
                except ValueError:
                    rate += eval(y[tag][p].replace("@0", "125")) * float(self.luminosity) / lumi

        return rate

    def __float__(self):
        return self.getrate()

__cache = {}
def getrate(*args):
    rate = __Rate(*args)
    if rate not in __cache:
        __cache[rate] = float(rate)
    return __cache[rate]

def getrates(*args, **kwargs):
    fmt = "rate {} {} {} {}"
    for kw, kwarg in kwargs.iteritems():
        if kw == "format":
            fmt = kwarg
    ggH, qqZZ, ggZZ, VBFbkg, ZX = getrate("ggH", *args), getrate("qqZZ", *args), getrate("ggZZ", *args), getrate("VBF bkg", *args), getrate("ZX", *args)

    result =  fmt.format(ggH, qqZZ, ggZZ, ZX)
    return result

def gettemplate(*args):
    try:
        try:
            return Template("prod+dec", *args).gettemplate()
        except ValueError as e1:
            return IntTemplate("prod+dec", *args).gettemplate()
    except ValueError as e2:
        raise ValueError("Can't gettemplate using args:\n{}\n\nTrying to make a regular template:\n{}\n\nTrying to make an interference template:\n{}".format(args, e1.message, e2.message))

def getdatatree(*args):
    return tfiles[DataTree(*args).treefile].candTree
    #it's empty if we don't actually unblind

def getsubtractdatatree(*args):
    return tfiles[SubtractDataTree(*args).treefile].candTree

def discriminantnames(*args):
    theset = set()
    for production in config.productionsforcombine:
        result = tuple(d.name for d in Template("ggH", "0+", "2e2mu", "prod+dec", production, *args).discriminants)
        theset.add(result)
    assert len(theset) == 1  #sanity check
    return result

class SigmaIOverSigma1(MultiEnum):
    enums = (Analysis, ProductionMode)
    def check(self, *args):
        if self.productionmode == "ggH":
            raise ValueError("need separate implementation for ggH, gi is defined differently for the pure sample\n{}".format(args))
    @property
    def result(self):
        sigmai = ReweightingSample(self.analysis.purehypotheses[1], self.productionmode).xsec
        sigma1 = ReweightingSample("SM", self.productionmode).xsec
        return sigmai/sigma1

def sigmaioversigma1(*args):
    return SigmaIOverSigma1(*args).result
