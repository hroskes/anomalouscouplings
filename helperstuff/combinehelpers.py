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
        if self.productionmode != "VBF bkg":
            return self.yamlrate()

        return Template("fa3", "prod+dec", "D_int_decay", self.category, self.productionmode, self.channel, self.production).gettemplate().Integral()*float(self.luminosity)

    def yamlrate(self):
        if self.production.year == 2016:
            filename = os.path.join(config.repositorydir, "helperstuff", "Datacards13TeV_ICHEP2016", "LegoCards", "configs", "inputs", "yields_per_tag_category_13TeV_{}.yaml".format(self.channel))
            tags = [tag.replace("Ichep16", "") for tag in self.category.names if "Ichep16" in tag]
        elif self.production.year == 2015:
            assert False
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
            productionmodes = ["ggH", "ttH"]
        elif self.productionmode == "VBF":
            productionmodes = ["qqH"]
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
    fmt = "rate {} {} {} {} {} {} {} {}"
    disableproductionmodes = ()
    useproductionmodes = ("ggH", "VBF", "WH", "ZH", "qqZZ", "ggZZ", "VBF bkg", "ZX")
    for kw, kwarg in kwargs.iteritems():
        if kw == "format":
            fmt = kwarg
        elif kw == "productionmodes":
            useproductionmodes = kwarg
        elif kw == "disableproductionmodes":
            disableproductionmodes = kwarg
        else:
            raise ValueError("Unknown kwarg {}={}!".format(kw, kwarg))

    rates = [getrate(c, *args) if c not in disableproductionmodes else 0 for c in useproductionmodes]

    result =  fmt.format(*rates)
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
    @property
    def result(self):
        if self.productionmode == "ggH":
            #for ggH gi for the pure BSM sample is defined so the xsecs are equal
            gi = getattr(
                         ReweightingSample(self.analysis.purehypotheses[1], self.productionmode),
                         self.analysis.couplingname
                        )
            return 1/gi**2
        if self.productionmode in ("VBF", "ZH", "WH"):
            #for VBF, ZH, and WH, gi for the pure BSM sample is defined = 1
            sigmai = ReweightingSample(self.analysis.purehypotheses[1], self.productionmode).xsec
            sigma1 = ReweightingSample("SM", self.productionmode).xsec
            return sigmai/sigma1
        assert False

def sigmaioversigma1(*args):
    return SigmaIOverSigma1(*args).result
