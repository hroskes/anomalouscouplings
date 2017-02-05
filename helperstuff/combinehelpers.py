import collections
import config
from enums import Analysis, categories, Category, Channel, EnumItem, MultiEnum, MyEnum, Production, ProductionMode
from math import copysign, isnan
import os
import ROOT
from samples import ReweightingSample, Sample
from templates import DataTree, IntTemplate, Template, TemplatesFile
from utilities import tfiles
import yaml

datacardprocessline = "ggH qqH WH ZH bkg_qqzz bkg_ggzz bkg_vbf bkg_zjets"
datacardprocessorder = [ProductionMode(p) for p in datacardprocessline.split()]

class LuminosityType(MyEnum):
    enumname = "luminositytype"
    enumitems = (
                 EnumItem("fordata"),
                 EnumItem("forexpectedscan"),
                 EnumItem("customluminosity"),
                )
    def __init__(self, value):
        tryfloat = True
        try:
            Production(value)
            tryfloat = False
        except ValueError:
            pass

        if (isinstance(value, MyEnum) or isinstance(value, MultiEnum)) and not (
              isinstance(value, LuminosityType) and value == "customluminosity"
           ):
            tryfloat = False

        if tryfloat:
            try:
                value = float(value)
                isfloat = True
            except (TypeError, ValueError):
                isfloat = False
        if tryfloat and isfloat:
            super(LuminosityType, self).__init__("customluminosity")
            self.customluminosity = value
        else:
            super(LuminosityType, self).__init__(value)

    def __float__(self):
        try:
            return self.customluminosity
        except AttributeError:
            raise TypeError("Can't convert LuminosityType('{}') to float".format(self))

    def __eq__(self, other):
        if isinstance(other, LuminosityType) and self == "customluminosity" == other:
            return float(other) == float(self)
        return super(LuminosityType, self).__eq__(other)

class Luminosity(MultiEnum):
    enumname = "luminosity"
    enums = [LuminosityType, Production]
    def check(self, *args):
        dontcheck = []
        if self.luminositytype != "fordata":
            dontcheck.append(Production)
        super(Luminosity, self).check(self, *args, dontcheck=dontcheck)
    def __float__(self):
        if self.luminositytype == "forexpectedscan": return float(config.expectedscanluminosity)
        if self.luminositytype == "fordata": return float(self.production.dataluminosity)
        if self.luminositytype == "customluminosity": return float(self.luminositytype)

class __Rate(MultiEnum):
    enums = [ProductionMode, Channel, Luminosity, Production, Category, Analysis]

    def getrate(self):
        if self.productionmode == "ZX" and not config.usedata:
            return self.yamlrate()

        if self.productionmode.issignal:
            hypothesis = self.analysis.purehypotheses[0]
        else:
            hypothesis = None
        return gettemplate(self.analysis, self.category, self.productionmode, self.channel, self.production, hypothesis).Integral()*float(self.luminosity)

    def yamlrate(self):
        if self.production.year == 2016:
            filename = os.path.join(config.repositorydir, "helperstuff", "Datacards13TeV_ICHEP2016", "LegoCards", "configs", "inputs", "yields_per_tag_category_13TeV_{}.yaml".format(self.channel))
            tags = [tag.replace("Ichep16", "").replace("tagged", "Tagged") for tag in self.category.names if "Ichep16" in tag]
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
        rate = 0

        for tag in tags:
            for p in self.productionmode.yamlratenames:
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
    disableproductionmodes = ()
    useproductionmodes = datacardprocessorder
    for kw, kwarg in kwargs.iteritems():
        if kw == "format":
            fmt = kwarg
        elif kw == "productionmodes":
            useproductionmodes = kwarg
        elif kw == "disableproductionmodes":
            disableproductionmodes = list(kwarg)
        else:
            raise ValueError("Unknown kwarg {}={}!".format(kw, kwarg))

    fmt = "rate " + " ".join("{}" for _ in useproductionmodes)

    rates = [getrate(c, *args) if c not in disableproductionmodes else 0 for c in useproductionmodes]

    result =  fmt.format(*rates)
    return result

def gettemplate(*args):
    try:
        try:
            t = Template(*args)
        except ValueError as e1:
            t = IntTemplate(*args)
    except ValueError as e2:
        raise ValueError("Can't gettemplate using args:\n{}\n\nTrying to make a regular template:\n{}\n\nTrying to make an interference template:\n{}".format(args, e1.message, e2.message))

    #############
    from datetime import date
    if date.today() > date(2017, 2, 8):
        raise ValueError("VBF bkg 4mu??!")
    if t.channel == "4mu" and t.productionmode == "VBFbkg":
        return gettemplate(t.analysis, t.category, t.productionmode, t.production, "4e")
    #############

    result = t.gettemplate()
    if isnan(result.Integral()):
        raise ValueError("{!r} has integral of nan!".format(t))

    return result

def getdatatree(*args):
    return tfiles[DataTree(*args).treefile].candTree
    #it's empty if we don't unblind

def getsubtractdatatree(*args):
    return tfiles[SubtractDataTree(*args).treefile].candTree

def getnobserved(*args):
    return getdatatree(*args).GetEntries()

def discriminants(*args):
    theset = set()
    for production in config.productionsforcombine:
        templatesfile = TemplatesFile("ggh", "2e2mu", production, *args)
        SM = templatesfile.analysis.purehypotheses[0]
        result = tuple(d for d in Template("ggH", SM, "2e2mu", production, *args).discriminants)
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
            sigma1 = ReweightingSample(self.analysis.purehypotheses[0], self.productionmode).xsec
            return sigmai/sigma1
        assert False

class SigmaIOverSigma1_VH(MultiEnum):
    enums = (Analysis,)
    @property
    def result(self):
        sigmai = sum(ReweightingSample(self.analysis.purehypotheses[1], VH).xsec for VH in ("ZH", "WH"))
        sigma1 = sum(ReweightingSample(self.analysis.purehypotheses[0], VH).xsec for VH in ("ZH", "WH"))
        return sigmai/sigma1

def sigmaioversigma1(*args):
    if "VH" in [str(_) for _ in args]:
        args = list(args)
        args = [_ for _ in args if str(_) != "VH"]
        return SigmaIOverSigma1_VH(*args).result
    return SigmaIOverSigma1(*args).result

def mixturesign(analysis, productionmode=None):
    analysis = Analysis(analysis)
    if productionmode is None:
        assert mixturesign(analysis, "ggH") == mixturesign(analysis, "VBF") == mixturesign(analysis, "ZH")
        productionmode = "ggH"
    return copysign(1, getattr(ReweightingSample(productionmode, analysis.mixdecayhypothesis), analysis.couplingname))
