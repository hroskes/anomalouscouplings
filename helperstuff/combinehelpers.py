import collections
from math import copysign, isnan
import os
import yaml

import ROOT

import config
from enums import Analysis, categories, Category, Channel, EnumItem, MultiEnum, MyEnum, Production, ProductionMode, shapesystematics
from samples import ReweightingSample, Sample
from templates import DataTree, IntTemplate, Template, TemplatesFile
from utilities import cache, tfiles

datacardprocessline = "ggH qqH WH ZH bkg_qqzz bkg_ggzz bkg_vbf bkg_zjets"
datacardprocessorder = [ProductionMode(p) for p in datacardprocessline.split()]

class LuminosityType(MyEnum):
    enumname = "luminositytype"
    enumitems = (
                 EnumItem("fordata"),
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
        if self.luminositytype == "fordata": return float(self.production.dataluminosity)
        if self.luminositytype == "customluminosity": return float(self.luminositytype)

class __Rate(MultiEnum):
    enums = [ProductionMode, Channel, Luminosity, Production, Category, Analysis]

    @cache
    def getrate(self):
        if self.copyfromotherrate: return self.copyfromotherrate.getrate()
        from yields import YieldValue
        result = YieldValue(self.production, self.productionmode, self.channel, self.category, self.analysis).value
        if self.productionmode == "ZX":
            if self.luminositytype != "fordata":
                result *= float(self.luminosity) / float(Luminosity("fordata", self.production))
        else:
            result *= float(self.luminosity)
        return result

    def __float__(self):
        return self.getrate()

    @property
    def copyfromotherrate(self):
        return None

def getrate(*args):
    return float(__Rate(*args))

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

    if t.shapesystematic not in shapesystematics:
        kwargs = {enum.enumname: getattr(t, enum.enumname) for enum in type(t).needenums}
        kwargs["shapesystematic"] = ""
        t = type(t)(*kwargs.values())  #can't use actual kwargs

    result = t.gettemplate()
    if isnan(result.Integral()):
        raise ValueError("{!r} has integral of nan!".format(t))

    return result

@cache
def zerotemplate(*args):
    result = gettemplate(*args)
    result = result.Clone("zerotemplate"+str(hash(args)))
    result.Reset("MICES")
    assert result.Integral() == 0
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
    templatesfile = TemplatesFile("ggh", "2e2mu", *args)
    SM = templatesfile.analysis.purehypotheses[0]
    result = tuple(d for d in Template("ggH", SM, "2e2mu", *args).discriminants)
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
            sigmai = ReweightingSample(self.analysis.purehypotheses[1], self.productionmode).nominalJHUxsec
            sigma1 = ReweightingSample(self.analysis.purehypotheses[0], self.productionmode).nominalJHUxsec
            return sigmai/sigma1
        assert False

class SigmaIOverSigma1_VH(MultiEnum):
    enums = (Analysis,)
    @property
    def result(self):
        import constants
        if self.analysis == "fa3":
            coupling = constants.g4VH
        elif self.analysis == "fa2":
            coupling = constants.g2VH
        elif self.analysis == "fL1":
            coupling = constants.g1prime2VH
        elif self.analysis == "fL1Zg":
            coupling = constants.ghzgs1prime2VH_gen
        else:
            assert False, self.analysis
        return 1/coupling**2 * sigmaioversigma1("ggH", self.analysis)

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
