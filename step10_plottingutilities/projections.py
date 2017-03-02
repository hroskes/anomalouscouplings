#!/usr/bin/env python
import abc
import collections
import glob
import itertools
from math import copysign, sqrt
import os
import re
import subprocess
from tempfile import mkdtemp

import ROOT

from helperstuff import config, constants, run1info, stylefunctions as style
from helperstuff.combinehelpers import getrate, gettemplate, Luminosity
from helperstuff.enums import analyses, Analysis, categories, Category, Channel, channels, EnumItem, MultiEnum, MultiEnumABCMeta, MyEnum, Production, ProductionMode, productions, ShapeSystematic
from helperstuff.samples import ReweightingSample, samplewithfai
import helperstuff.rootoverloads.histogramaxisnumbers, helperstuff.rootoverloads.histogramfloor
from helperstuff.templates import IntTemplate, Template, TemplatesFile
from helperstuff.utilities import cache, tfiles, pairwise

c1 = ROOT.TCanvas("c1", "", 8, 30, 800, 800)
style.applycanvasstyle(c1)

class TemplateForProjection(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, *args, **kwargs):
        self.linecolor = 0
        self.linestyle = 1
        self.linewidth = 1
        self.fillcolor = 0
        self.fillstyle = 0
        self.legendoption = "l"
        for kw, kwarg in kwargs.items():
          if kw in ("linecolor", "linestyle", "linewidth", "fillcolor", "fillstyle", "legendoption"):
            setattr(self, kw, kwarg)
            del kwargs[kw]
        self.__initedhistogram = False
        try:
            super(TemplateForProjection, self).__init__(*args, **kwargs)
        except TypeError:
            print args
            print kwargs
            raise

    @abc.abstractmethod
    def inithistogram(self):
        self.__initedhistogram = True
        for attr in self.histogramattrs:
            assert hasattr(self, attr)

    @cache
    def Projection(self, i, rebin=None):
        if rebin is not None:
            result = self.Projection(i).Clone()
            result.Rebin(rebin)
            return result
        result = self.h.Projection(i)
        result.SetName(self.discriminants[i].name)
        result.SetTitle(self.discriminants[i].title)
        result.SetXTitle(self.discriminants[i].title)
        result.SetLineColor(self.linecolor)
        result.SetLineStyle(self.linestyle)
        result.SetLineWidth(self.linewidth)
        result.SetFillColor(self.fillcolor)
        result.SetFillStyle(self.fillstyle)
        return result

    def AddToLegend(self, legend):
        legend.AddEntry(self.Projection(0), self.title, self.legendoption)

    def Integral(self):
        if self.cantakeintegral:
            return self.h.Integral()
        raise ValueError("Can't take integral of {}".format(self))

    def Floor(self, *args, **kwargs):
        self.h.Floor(*args, **kwargs)

    @abc.abstractproperty
    def discriminants(self):
        pass

    def __getattr__(self, attr):
        if not self.__initedhistogram and attr in self.histogramattrs:
            self.inithistogram()
            return getattr(self, attr)
        raise AttributeError("'{}' object has no attribute '{}'".format(type(self).__name__, attr))

    @property
    def histogramattrs(self):
        return "h", "hstackoption", "cantakeintegral"

class Normalization(MyEnum):
    enumname = "normalization"
    enumitems = [
                 EnumItem("rescalemixtures"),
                ]
normalizations = Normalization.items()

class EnrichStatus(MyEnum):
    enumname = "enrichstatus"
    enumitems = [
                 EnumItem("enrich"),
                 EnumItem("fullrange", "noenrich"),
                 EnumItem("blind", "impoverish"),
                ]
    def cuttext(self):
        if self == "enrich": return "D_{bkg} > 0.5"
        if self == "noenrich": return ""
        if self == "impoverish": return "D_{bkg} < 0.5"
        assert False
    def dirname(self):
        return str(self)
enrichstatuses = EnrichStatus.items()

class BaseTemplateFromFile(TemplateForProjection):
    def __init__(self, *args, **kwargs):
        super(BaseTemplateFromFile, self).__init__(*args, **kwargs)
        self.didenrich = False
        #I don't know what the following line accomplishes
        # but it fixes the WH interference templates.
        #Otherwise they are wrong in some way that
        # I don't understand
        self.template.gettemplate()

    def inithistogram(self):
        self.h = self.template.gettemplate().Clone()

        if self.productionmode in ["ggH", "VBF", "ZH", "WH"]:
            scalefactor = getrate(self.channel, self.category, self.productionmode, "fordata", self.production, self.analysis) / Template(self.production, self.category, self.analysis, self.channel, self.productionmode, self.analysis.purehypotheses[0]).gettemplate().Integral()
        elif self.productionmode == "data" or self.h.Integral() == 0:
            scalefactor = 1
        else:
            scalefactor = getrate(self.channel, self.category, self.productionmode, "fordata", self.production, self.analysis) / self.h.Integral()

        if self.productionmode == "data":
            self.hstackoption = "P"
            self.h.SetMarkerColor(self.linecolor)
            self.h.SetMarkerStyle(20)
        else:
            self.hstackoption = "hist"

        self.h.Scale(scalefactor)

        self.doenrich()

        self.cantakeintegral = (self.enrichstatus == "fullrange")

        super(BaseTemplateFromFile, self).inithistogram()

    def check(self, *args):
        if self.normalization is None:
            self.normalization = Normalization("rescalemixtures")
        super(BaseTemplateFromFile, self).check(*args)

    @abc.abstractproperty
    def title(self):
        pass

    @property
    def discriminants(self):
        return self.templatesfile.discriminants

    def doenrich(self):
        if self.didenrich: return
        self.didenrich = True
        for x, y, z in itertools.product(range(1, self.h.GetNbinsX()+1), range(1, self.h.GetNbinsY()+1), range(1, self.h.GetNbinsZ()+1)):
            if (self.enrichstatus == "impoverish" and self.h.GetZaxis().GetBinLowEdge(z) >= .5 \
             or self.enrichstatus == "enrich"     and self.h.GetZaxis().GetBinLowEdge(z) < .5):
                self.h.SetBinContent(x, y, z, 0)

class TemplateFromFile(BaseTemplateFromFile, MultiEnum):
    __metaclass__ = MultiEnumABCMeta
    enums = [Template, Normalization, EnrichStatus]
    @property
    def title(self):
        return self.template.title()

class IntTemplateFromFile(BaseTemplateFromFile, MultiEnum):
    __metaclass__ = MultiEnumABCMeta
    enums = [IntTemplate, Normalization, EnrichStatus]
    @property
    def template(self):
        return self.inttemplate
    @property
    def title(self):
        return ""

class TemplateSumBase(TemplateForProjection):
    def __init__(self, title, *templatesandfactors, **kwargs):
        self.title = title
        self.templatesandfactors = templatesandfactors
        analyses = {t[0].analysis for t in templatesandfactors}
        assert len(analyses) == 1
        self.analysis = analyses.pop()
        super(TemplateSumBase, self).__init__(**kwargs)

    def __hash__(self):
        return hash((self.title, self.templatesandfactors))

    def __eq__(self, other):
        return hash(self) == hash(other)

    @abc.abstractmethod
    def inithistogram(self):
        super(TemplateSumBase, self).inithistogram()

    @property
    def discriminants(self):
        assert len({template.discriminants for template, factor in self.templatesandfactors}) == 1
        return self.templatesandfactors[0][0].discriminants

class TemplateSum(TemplateSumBase):
    def inithistogram(self):
        self.cantakeintegral = True
        self.h = None
        for template, factor in self.templatesandfactors:
            if not factor: continue
            if self.h is None:
                self.h = template.h.Clone()
                self.h.Scale(factor)
            else:
                self.h.Add(template.h, factor)
            self.cantakeintegral = self.cantakeintegral and template.cantakeintegral
        self.hstackoption = "hist"
        super(TemplateSum, self).inithistogram()

class IntegralSum(TemplateSumBase):
    counter = 0
    def __init__(self, *templatesandfactors, **kwargs):
        if "linecolor" not in kwargs: kwargs.update(linecolor=0)
        super(IntegralSum, self).__init__("", *templatesandfactors, **kwargs)
    @classmethod
    def hname(cls):
        cls.counter += 1
        return "h{}".format(cls.counter)
    def inithistogram(self):
        self.hstackoption = None
        integral = 0
        self.cantakeintegral = True
        for template, factor in self.templatesandfactors:
            if not factor: continue
            self.cantakeintegral = self.cantakeintegral and template.cantakeintegral
            integral += template.Integral() * factor
        self.h = ROOT.TH1D(self.hname(), "", 1, 0, 1)
        self.h.Fill(.5, integral)
        assert abs(self.h.Integral() - integral)/integral < 1e-10, (self.h.Integral(), integral)
        super(IntegralSum, self).inithistogram()

class ComponentTemplateSum(TemplateSum):
    def __init__(self, title, SMintegral, *templatesandfactors, **kwargs):
        """
        Works the same as TemplateSum, but rescales
        """
        self.SMintegral = SMintegral
        super(ComponentTemplateSum, self).__init__(title, *templatesandfactors, **kwargs)
        for template, factor in templatesandfactors:
            if not isinstance(template, BaseTemplateFromFile):
                raise TypeError("{} can only come from TemplatesFromFiles".format(type(self).__name__))

    def __hash__(self):
        return hash((self.SMintegral, super(ComponentTemplateSum, self).__hash__()))

    def inithistogram(self):
        super(ComponentTemplateSum, self).inithistogram()
        if isinstance(self.SMintegral, TemplateForProjection):
            if not self.cantakeintegral:
                raise ValueError("Have to be able to take integral of ComponentTemplateSum, or normalize to a number")
            self.SMintegral = self.SMintegral.Integral()
        if self.Integral():
            self.h.Scale(self.SMintegral / self.h.Integral())

class ComponentTemplateSumInGroup(TemplateSum):
    def __init__(self, title, mygroup, SMgroup, *templatesandfactors, **kwargs):
        self.mygroup = mygroup
        self.SMgroup = SMgroup
        super(ComponentTemplateSumInGroup, self).__init__(title, *templatesandfactors, **kwargs)
        for template, factor in templatesandfactors:
            if not isinstance(template, BaseTemplateFromFile):
                raise TypeError("{} can only come from TemplatesFromFiles".format(type(self).__name__))

    def __hash__(self):
        return hash((self.mygroup, self.SMgroup, super(ComponentTemplateSumInGroup, self).__hash__()))

    def inithistogram(self):
        super(ComponentTemplateSumInGroup, self).inithistogram()
        self.h.Scale(self.SMgroup.Integral() / self.mygroup.Integral())

class Projections(MultiEnum):
  enums = [Analysis, Normalization, ShapeSystematic, Production, EnrichStatus]
  enumname = "projections"

  def applysynonyms(self, enumsdict):
    super(Projections, self).applysynonyms(enumsdict)

  def check(self, *args):
    dontcheck = []
    if self.normalization is None:
      self.normalization = Normalization("rescalemixtures")
    if self.shapesystematic is None:
      self.shapesystematic = ShapeSystematic("")
    super(Projections, self).check(*args, dontcheck=dontcheck)

  @cache
  def TemplateFromFile(self, *args, **kwargs):
    return TemplateFromFile(*args, **kwargs)
  @cache
  def IntTemplateFromFile(self, *args, **kwargs):
    return IntTemplateFromFile(*args, **kwargs)
  @cache
  def TemplateSum(self, *args, **kwargs):
    return TemplateSum(*args, **kwargs)
  @cache
  def IntegralSum(self, *args, **kwargs):
    return IntegralSum(*args, **kwargs)
  @cache
  def ComponentTemplateSum(self, *args, **kwargs):
    return ComponentTemplateSum(*args, **kwargs)
  @cache
  def ComponentTemplateSumInGroup(self, *args, **kwargs):
    return ComponentTemplateSumInGroup(*args, **kwargs)

  class CategoryAndChannel(MultiEnum): enums = (Category, Channel)

  def projections(self, *categoryandchannel, **kwargs):
    info = self.CategoryAndChannel(*categoryandchannel)
    category, channel = info.category, info.channel

    giname = self.analysis.couplingname
    BSMhypothesis = self.analysis.purehypotheses[1]

    ggHfactor = VBFfactor = VHfactor = ttHfactor = 1
    subdir = ""
    saveasdir = self.saveasdir(info)
    customfai = None
    animation = False
    saveasappend = ""
    exts = "png", "eps", "root", "pdf"
    otherthingstodraw = []
    otherthingstoaddtolegend = []
    floor = False
    nicestyle = False
    justoneproductionmode = None
    legendargs = [(.65, .6, .9, .9)]*3
    rebin = None
    muV = muf = 1
    for kw, kwarg in kwargs.iteritems():
       if kw == "productionmode":
           justoneproductionmode = True
           assert "muVmuf" not in kwargs
           if   kwarg == "ggH": ggHfactor = 1; VBFfactor = VHfactor  = ttHfactor = 0
           elif kwarg == "VBF": VBFfactor = 1; ggHfactor = VHfactor  = ttHfactor = 0
           elif kwarg ==  "VH": VHfactor  = 1; ggHfactor = VBFfactor = ttHfactor = 0
           elif kwarg == "ttH": ttHfactor = 1; ggHfactor = VVFfactor = VHfactor  = 0
           else: raise ValueError("Unknown productionmode {}!".format(kwarg))
       elif kw == "muVmuf":
           assert "productionmode" not in kwargs
           muV, muf = kwarg
       elif kw == "subdir":
           subdir = kwarg
       elif kw == "saveasdir":
           saveasdir = kwarg
       elif kw == "saveasappend":
           saveasappend = kwarg
       elif kw == "customfaiforanimation":
           animation = True
           customfai, customfaiproductionmode = kwarg
           customfai = float(customfai)
           customfaiproductionmode = ProductionMode(customfaiproductionmode)
           customsample = samplewithfai(customfaiproductionmode, self.analysis, customfai)
           g1_custom = customsample.g1
           gi_custom = getattr(customsample, BSMhypothesis.couplingname)
           exts = "gif",
       elif kw == "otherthingstodraw":
           otherthingstodraw += kwarg
       elif kw == "floor":
           floor = bool(int(kwarg))
       elif kw == "nicestyle":
           nicestyle = bool(int(kwarg))
           assert saveasdir not in kwargs
           saveasdir = self.saveasdir_niceplots(category)
           if nicestyle and channel != "2e2mu": return
       else:
           raise TypeError("Unknown kwarg {}={}".format(kw, kwarg))

    SMhypothesis = self.analysis.purehypotheses[0]
    gi_ggHBSM = getattr(ReweightingSample("ggH", BSMhypothesis), BSMhypothesis.couplingname)
    gi_VBFBSM = copysign((ReweightingSample("VBF", SMhypothesis).xsec / ReweightingSample("VBF", BSMhypothesis).xsec)**.25, gi_ggHBSM)
    gi_VHBSM = copysign(((ReweightingSample("WH", SMhypothesis).xsec + ReweightingSample("ZH", SMhypothesis).xsec) / (ReweightingSample("WH", BSMhypothesis).xsec + ReweightingSample("ZH", BSMhypothesis).xsec))**.25, gi_ggHBSM)
    gi_VVHBSM = copysign(((ReweightingSample("VBF", SMhypothesis).xsec + ReweightingSample("WH", SMhypothesis).xsec + ReweightingSample("ZH", SMhypothesis).xsec) / (ReweightingSample("VBF", BSMhypothesis).xsec + ReweightingSample("WH", BSMhypothesis).xsec + ReweightingSample("ZH", BSMhypothesis).xsec))**.25, gi_ggHBSM)
    if category == "UntaggedMor17":
        g1_mix = 1/sqrt(2)
        gi_mix = 1/sqrt(2)*gi_ggHBSM
        fainame = self.analysis.title(superscript="dec")
    elif category == "VBF2jTaggedMor17":
        g1_mix = 1/2**.25
        gi_mix = 1/2**.25 * gi_VBFBSM
        fainame = self.analysis.title(superscript="VBFdec")
    elif category == "VHHadrTaggedMor17":
        g1_mix = 1/2**.25
        gi_mix = 1/2**.25 * gi_VHBSM
        fainame = self.analysis.title(superscript="VHdec")

    allggHg12gi0 = {}
    allggHg10gi2 = {}
    allggHg11gi1 = {}
    allttHg12gi0 = {}
    allttHg10gi2 = {}
    allttHg11gi1 = {}

    allVBFg14gi0 = {}
    allVBFg13gi1 = {}
    allVBFg12gi2 = {}
    allVBFg11gi3 = {}
    allVBFg10gi4 = {}
    allVBFpieces = {}

    allZHg14gi0 = {}
    allZHg13gi1 = {}
    allZHg12gi2 = {}
    allZHg11gi3 = {}
    allZHg10gi4 = {}
    allZHpieces = {}

    allWHg14gi0 = {}
    allWHg13gi1 = {}
    allWHg12gi2 = {}
    allWHg11gi3 = {}
    allWHg10gi4 = {}
    allWHpieces = {}

    for ca, ch in itertools.product(categories, channels):
      allggHg12gi0[ca,ch] = self.TemplateFromFile(   "ggH", ca, "fullrange", self.normalization, self.production, ch, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
      allggHg10gi2[ca,ch] = self.TemplateFromFile(   "ggH", ca, "fullrange", self.normalization, self.production, ch, self.shapesystematic, self.analysis, BSMhypothesis)
      allggHg11gi1[ca,ch] = self.IntTemplateFromFile("ggH", ca, "fullrange", self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g11gi1")

      allttHg12gi0[ca,ch] = self.TemplateFromFile(   "ttH", ca, "fullrange", self.normalization, self.production, ch, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
      allttHg10gi2[ca,ch] = self.TemplateFromFile(   "ttH", ca, "fullrange", self.normalization, self.production, ch, self.shapesystematic, self.analysis, BSMhypothesis)
      allttHg11gi1[ca,ch] = self.IntTemplateFromFile("ttH", ca, "fullrange", self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g11gi1")

      allVBFg14gi0[ca,ch] = self.TemplateFromFile(   "VBF", ca, "fullrange", self.normalization, self.production, ch, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
      allVBFg13gi1[ca,ch] = self.IntTemplateFromFile("VBF", ca, "fullrange", self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g13gi1")
      allVBFg12gi2[ca,ch] = self.IntTemplateFromFile("VBF", ca, "fullrange", self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g12gi2")
      allVBFg11gi3[ca,ch] = self.IntTemplateFromFile("VBF", ca, "fullrange", self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g11gi3")
      allVBFg10gi4[ca,ch] = self.TemplateFromFile(   "VBF", ca, "fullrange", self.normalization, self.production, ch, self.shapesystematic, self.analysis, BSMhypothesis)

      allVBFpieces[ca,ch] = [allVBFg14gi0[ca,ch], allVBFg13gi1[ca,ch], allVBFg12gi2[ca,ch], allVBFg11gi3[ca,ch], allVBFg10gi4[ca,ch]]

      allZHg14gi0[ca,ch] = self.TemplateFromFile(   "ZH", ca, "fullrange", self.normalization, self.production, ch, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
      allZHg13gi1[ca,ch] = self.IntTemplateFromFile("ZH", ca, "fullrange", self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g13gi1")
      allZHg12gi2[ca,ch] = self.IntTemplateFromFile("ZH", ca, "fullrange", self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g12gi2")
      allZHg11gi3[ca,ch] = self.IntTemplateFromFile("ZH", ca, "fullrange", self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g11gi3")
      allZHg10gi4[ca,ch] = self.TemplateFromFile(   "ZH", ca, "fullrange", self.normalization, self.production, ch, self.shapesystematic, self.analysis, BSMhypothesis)

      allZHpieces[ca,ch] = [allZHg14gi0[ca,ch], allZHg13gi1[ca,ch], allZHg12gi2[ca,ch], allZHg11gi3[ca,ch], allZHg10gi4[ca,ch]]

      allWHg14gi0[ca,ch] = self.TemplateFromFile(   "WH", ca, "fullrange", self.normalization, self.production, ch, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
      allWHg13gi1[ca,ch] = self.IntTemplateFromFile("WH", ca, "fullrange", self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g13gi1")
      allWHg12gi2[ca,ch] = self.IntTemplateFromFile("WH", ca, "fullrange", self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g12gi2")
      allWHg11gi3[ca,ch] = self.IntTemplateFromFile("WH", ca, "fullrange", self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g11gi3")
      allWHg10gi4[ca,ch] = self.TemplateFromFile(   "WH", ca, "fullrange", self.normalization, self.production, ch, self.shapesystematic, self.analysis, BSMhypothesis)

      allWHpieces[ca,ch] = [allWHg14gi0[ca,ch], allWHg13gi1[ca,ch], allWHg12gi2[ca,ch], allWHg11gi3[ca,ch], allWHg10gi4[ca,ch]]

    ffHSM     = self.IntegralSum(*sum(
                                      ([(allggHg12gi0[ca,ch], 1), (allttHg12gi0[ca,ch], 1)]
                                          for ca, ch in itertools.product(categories, channels)),
                                       []
                                      ))
    ffHBSM    = self.IntegralSum(*sum(
                                      ([(allggHg10gi2[ca,ch], 1), (allttHg10gi2[ca,ch], 1)]
                                         for ca, ch in itertools.product(categories, channels)),
                                      []
                                     ))
    ffHmix_p  = self.IntegralSum(*sum(
                                      ([(allggHg12gi0[ca,ch], g1_mix**2), (allggHg10gi2[ca,ch], (gi_mix/gi_ggHBSM)**2), (allggHg11gi1[ca,ch],  g1_mix*gi_mix/gi_ggHBSM),
                                        (allttHg12gi0[ca,ch], g1_mix**2), (allttHg10gi2[ca,ch], (gi_mix/gi_ggHBSM)**2), (allttHg11gi1[ca,ch],  g1_mix*gi_mix/gi_ggHBSM)]
                                         for ca, ch in itertools.product(categories, channels)),
                                       []
                                      ))
    ffHmix_m  = self.IntegralSum(*sum(
                                      ([(allggHg12gi0[ca,ch], g1_mix**2), (allggHg10gi2[ca,ch], (gi_mix/gi_ggHBSM)**2), (allggHg11gi1[ca,ch], -g1_mix*gi_mix/gi_ggHBSM),
                                        (allttHg12gi0[ca,ch], g1_mix**2), (allttHg10gi2[ca,ch], (gi_mix/gi_ggHBSM)**2), (allttHg11gi1[ca,ch], -g1_mix*gi_mix/gi_ggHBSM)]
                                         for ca, ch in itertools.product(categories, channels)),
                                       []
                                     ))

    VVHSM     = self.IntegralSum(*sum(
                                      ([(allVBFg14gi0[ca,ch], 1), (allZHg14gi0[ca,ch], 1), (allWHg14gi0[ca,ch], 1)]
                                         for ca, ch in itertools.product(categories, channels)),
                                      []
                                     ))
    VVHBSM    = self.IntegralSum(*sum(
                                      ([(allVBFg10gi4[ca,ch], 1), (allZHg10gi4[ca,ch], 1), (allWHg10gi4[ca,ch], 1)]
                                         for ca, ch in itertools.product(categories, channels)),
                                      []
                                     ))
    VVHmix_p  = self.IntegralSum(*sum((
                                        [(template, g1_mix**(4-j) * (+gi_mix)**j) for j, template in enumerate(allVBFpieces[ca,ch])]
                                      + [(template, g1_mix**(4-j) * (+gi_mix)**j) for j, template in enumerate(allZHpieces[ca,ch])]
                                      + [(template, g1_mix**(4-j) * (+gi_mix)**j) for j, template in enumerate(allWHpieces[ca,ch])]
                                         for ca, ch in itertools.product(categories, channels)),
                                      []
                                     ))
    VVHmix_m  = self.IntegralSum(*sum((
                                        [(template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(allVBFpieces[ca,ch])]
                                      + [(template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(allZHpieces[ca,ch])]
                                      + [(template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(allWHpieces[ca,ch])]
                                         for ca, ch in itertools.product(categories, channels)),
                                      []
                                     ))

    del ca, ch

    ggHg12gi0 = {}
    ggHg10gi2 = {}
    ggHg11gi1 = {}
    ttHg12gi0 = {}
    ttHg10gi2 = {}
    ttHg11gi1 = {}

    VBFg14gi0 = {}
    VBFg13gi1 = {}
    VBFg12gi2 = {}
    VBFg11gi3 = {}
    VBFg10gi4 = {}
    VBFpieces = {}

    ZHg14gi0 = {}
    ZHg13gi1 = {}
    ZHg12gi2 = {}
    ZHg11gi3 = {}
    ZHg10gi4 = {}
    ZHpieces = {}

    WHg14gi0 = {}
    WHg13gi1 = {}
    WHg12gi2 = {}
    WHg11gi3 = {}
    WHg10gi4 = {}
    WHpieces = {}

    for ch in channels:
      ggHg12gi0[ch] = self.TemplateFromFile(   "ggH", category, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
      ggHg10gi2[ch] = self.TemplateFromFile(   "ggH", category, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, BSMhypothesis)
      ggHg11gi1[ch] = self.IntTemplateFromFile("ggH", category, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g11gi1")

      ttHg12gi0[ch] = self.TemplateFromFile(   "ttH", category, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
      ttHg10gi2[ch] = self.TemplateFromFile(   "ttH", category, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, BSMhypothesis)
      ttHg11gi1[ch] = self.IntTemplateFromFile("ttH", category, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g11gi1")

      VBFg14gi0[ch] = self.TemplateFromFile(   "VBF", category, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
      VBFg13gi1[ch] = self.IntTemplateFromFile("VBF", category, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g13gi1")
      VBFg12gi2[ch] = self.IntTemplateFromFile("VBF", category, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g12gi2")
      VBFg11gi3[ch] = self.IntTemplateFromFile("VBF", category, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g11gi3")
      VBFg10gi4[ch] = self.TemplateFromFile(   "VBF", category, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, BSMhypothesis)

      VBFpieces[ch] = [VBFg14gi0[ch], VBFg13gi1[ch], VBFg12gi2[ch], VBFg11gi3[ch], VBFg10gi4[ch]]

      ZHg14gi0[ch] = self.TemplateFromFile(   "ZH", category, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
      ZHg13gi1[ch] = self.IntTemplateFromFile("ZH", category, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g13gi1")
      ZHg12gi2[ch] = self.IntTemplateFromFile("ZH", category, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g12gi2")
      ZHg11gi3[ch] = self.IntTemplateFromFile("ZH", category, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g11gi3")
      ZHg10gi4[ch] = self.TemplateFromFile(   "ZH", category, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, BSMhypothesis)

      ZHpieces[ch] = [ZHg14gi0[ch], ZHg13gi1[ch], ZHg12gi2[ch], ZHg11gi3[ch], ZHg10gi4[ch]]

      WHg14gi0[ch] = self.TemplateFromFile(   "WH", category, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
      WHg13gi1[ch] = self.IntTemplateFromFile("WH", category, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g13gi1")
      WHg12gi2[ch] = self.IntTemplateFromFile("WH", category, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g12gi2")
      WHg11gi3[ch] = self.IntTemplateFromFile("WH", category, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g11gi3")
      WHg10gi4[ch] = self.TemplateFromFile(   "WH", category, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, BSMhypothesis)

      WHpieces[ch] = [WHg14gi0[ch], WHg13gi1[ch], WHg12gi2[ch], WHg11gi3[ch], WHg10gi4[ch]]

    del ch

    ggHSM = self.TemplateSum("ggH SM", (ggHg12gi0[channel], 1))
    ggHBSM = self.ComponentTemplateSumInGroup("ggH {}=1".format(fainame), ffHBSM, ffHSM, (ggHg10gi2[channel], 1))
    ggHmix_p  = self.ComponentTemplateSumInGroup("ggH {}=#plus0.5" .format(fainame), ffHmix_p, ffHSM,
                                                 (ggHg12gi0[channel], g1_mix**2), (ggHg10gi2[channel], (gi_mix/gi_ggHBSM)**2), (ggHg11gi1[channel],  g1_mix*gi_mix/gi_ggHBSM))
    ggHmix_m  = self.ComponentTemplateSumInGroup("ggH {}=#minus0.5".format(fainame), ffHmix_m, ffHSM,
                                                 (ggHg12gi0[channel], g1_mix**2), (ggHg10gi2[channel], (gi_mix/gi_ggHBSM)**2), (ggHg11gi1[channel], -g1_mix*gi_mix/gi_ggHBSM))
    ttHSM = self.TemplateSum("ttH SM", (ttHg12gi0[channel], 1))
    ttHBSM = self.ComponentTemplateSumInGroup("ttH {}=1".format(fainame), ffHBSM, ffHSM, (ttHg10gi2[channel], 1))
    ttHmix_p  = self.ComponentTemplateSumInGroup("ttH {}=#plus0.5" .format(fainame), ffHmix_p, ffHSM,
                                                 (ttHg12gi0[channel], g1_mix**2), (ttHg10gi2[channel], (gi_mix/gi_ggHBSM)**2), (ttHg11gi1[channel],  g1_mix*gi_mix/gi_ggHBSM))
    ttHmix_m  = self.ComponentTemplateSumInGroup("ttH {}=#minus0.5".format(fainame), ffHmix_m, ffHSM,
                                                 (ttHg12gi0[channel], g1_mix**2), (ttHg10gi2[channel], (gi_mix/gi_ggHBSM)**2), (ttHg11gi1[channel], -g1_mix*gi_mix/gi_ggHBSM))

    VBFSM = self.TemplateSum("VBF SM", (VBFg14gi0[channel], 1))
    VBFBSM = self.ComponentTemplateSumInGroup("VBF {}=1".format(fainame), VVHBSM, VVHSM, (VBFg10gi4[channel], 1))
    VBFmix_p  = self.ComponentTemplateSumInGroup("VBF {}=#plus0.5" .format(fainame), VVHmix_p, VVHSM,
                                            *((template, g1_mix**(4-j) * gi_mix**j) for j, template in enumerate(VBFpieces[channel]))
                                           )
    VBFmix_m  = self.ComponentTemplateSumInGroup("VBF {}=#minus0.5".format(fainame), VVHmix_m, VVHSM,
                                            *((template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(VBFpieces[channel]))
                                           )

    VHSM = self.TemplateSum("VH SM", (ZHg14gi0[channel], 1), (WHg14gi0[channel], 1))
    VHBSM = self.ComponentTemplateSumInGroup("VH {}=1".format(fainame), VVHBSM, VVHSM, (ZHg10gi4[channel], 1), (WHg10gi4[channel], 1))
    VHmix_p  = self.ComponentTemplateSumInGroup("VH {}=#plus0.5" .format(fainame), VVHmix_p, VVHSM,
                                            *(
                                                [(template, g1_mix**(4-j) * gi_mix**j) for j, template in enumerate(ZHpieces[channel])]
                                              + [(template, g1_mix**(4-j) * gi_mix**j) for j, template in enumerate(WHpieces[channel])]
                                             )
                                           )
    VHmix_m  = self.ComponentTemplateSumInGroup("VH {}=#minus0.5" .format(fainame), VVHmix_m, VVHSM,
                                            *(
                                                [(template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(ZHpieces[channel])]
                                              + [(template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(WHpieces[channel])]
                                             )
                                           )

    SM    = self.TemplateSum("SM",                                 (ggHSM,    ggHfactor), (VBFSM,    VBFfactor), (VHSM,    VHfactor), (ttHSM,    ttHfactor), linecolor=1)
    BSM   = self.TemplateSum("{}=1".format(self.analysis.title()), (ggHBSM,   ggHfactor), (VBFBSM,   VBFfactor), (VHBSM,   VHfactor), (ttHBSM,   ttHfactor), linecolor=ROOT.kCyan)
    mix_p = self.TemplateSum("{}=#plus0.5".format(fainame),        (ggHmix_p, ggHfactor), (VBFmix_p, VBFfactor), (VHmix_p, VHfactor), (ttHmix_p, ttHfactor), linecolor=ROOT.kGreen+3)
    mix_m = self.TemplateSum("{}=#minus0.5".format(fainame),       (ggHmix_m, ggHfactor), (VBFmix_m, VBFfactor), (VHmix_m, VHfactor), (ttHmix_m, ttHfactor), linecolor=4)

    qqZZ      = self.TemplateFromFile(category, self.enrichstatus, self.normalization, self.analysis, channel, "qqZZ",    self.shapesystematic, self.production, linecolor=6)
    ggZZ      = self.TemplateFromFile(category, self.enrichstatus, self.normalization, self.analysis, channel, "ggZZ",    self.shapesystematic, self.production, linecolor=ROOT.kOrange+6)
    VBFbkg    = self.TemplateFromFile(category, self.enrichstatus, self.normalization, self.analysis, channel, "VBF bkg", self.shapesystematic, self.production, linecolor=ROOT.kViolet+7)
    ZX        = self.TemplateFromFile(category, self.enrichstatus, self.normalization, self.analysis, channel, "ZX",      self.shapesystematic, self.production, linecolor=2)

    templates = []
    if not animation:
        templates += [
                      SM, BSM, mix_p, mix_m,
                     ]
    else:
        assert ggHfactor == VBFfactor == VHfactor == 1
        if customfai <  0: plusminus = "#minus"
        if customfai == 0: plusminus = ""
        if customfai  > 0: plusminus = "#plus"
        ggHcustom = self.ComponentTemplateSum("ggH ({}={}{:.2f})".format(self.analysis.title(superscript="dec"), plusminus, abs(customsample.fai("ggH", self.analysis))), 1, (ggHg12gi0[channel], g1_custom**2), (ggHg10gi2[channel], (gi_custom/gi_ggHBSM)**2), (ggHg11gi1[channel],  g1_custom*gi_custom/gi_ggHBSM), linecolor=1)
        VBFcustom = self.ComponentTemplateSum("VBF ({}={}{:.2f})".format(self.analysis.title(superscript="VBF"), plusminus, abs(customsample.fai("VBF", self.analysis))), 1,
                                              *((template, g1_custom**(4-j) * gi_custom**j) for j, template in enumerate(VBFpieces[channel])),
                                              linecolor=2
                                             )
        VHcustom  = self.ComponentTemplateSum("VH ({}={}{:.2f})".format(self.analysis.title(superscript="VH"), plusminus, abs(customsample.fai("VH", self.analysis))), 1,
                                              *(
                                                  [(template, g1_custom**(4-j) * gi_custom**j) for j, template in enumerate(ZHpieces[channel])]
                                                + [(template, g1_custom**(4-j) * gi_custom**j) for j, template in enumerate(WHpieces[channel])]
                                               ),
                                               linecolor=4
                                             )
        for t in ggHcustom, VBFcustom, VHcustom:
            assert abs(t.h.Integral() - 1) < 1e-6, t.h.Integral()

        templates += [
                      ggHcustom, VBFcustom, VHcustom
                     ]

    if not justoneproductionmode and not animation:
        templates += [
                      qqZZ, ggZZ, VBFbkg, ZX,
                     ]

    if nicestyle:
        ZX = self.TemplateSum("Z+X",
                              *((self.TemplateFromFile(category, self.enrichstatus, self.normalization, self.analysis, ch, "ZX",      self.shapesystematic, self.production, linecolor=2), 1) for ch in channels),
                              linecolor=1, fillcolor=ROOT.TColor.GetColor("#669966"), linewidth=2, fillstyle=1001, legendoption="f")
        ZZ = self.TemplateSum("ZZ/Z#gamma*",
                              *((self.TemplateFromFile(category, self.enrichstatus, self.normalization, self.analysis, ch, p,      self.shapesystematic, self.production, linecolor=2), 1) for ch in channels for p in ("qqZZ", "ggZZ", "VBFbkg", "ZX")),
                              linecolor=1, fillcolor=ROOT.kAzure-9, linewidth=2, fillstyle=1001, legendoption="f")
        if category == "Untagged":
            superscript = None
        elif category == "VBFtagged":
            superscript = "VBF"
        elif category == "VHHadrtagged":
            superscript = "VH"

        templates = [ZZ, ZX]

        legendargs = [(0.20,0.57,0.58,0.90), (0.20,0.57,0.58,0.90), (0.23,0.57,0.61,0.90)]
        if not animation:
            SMffH = self.TemplateSum("",
                                     *sum(
                                          ([(ggHg12gi0[ch], 1), (ttHg12gi0[ch], 1)]
                                              for ch in channels),
                                           []
                                         )
                                    )
            SMVVH = self.TemplateSum("",
                                     *sum(
                                          ([(VBFg14gi0[ch], 1), (ZHg14gi0[ch], 1), (WHg14gi0[ch], 1)]
                                              for ch in channels),
                                           []
                                         )
                                    )
            BSMffH = self.ComponentTemplateSumInGroup("", ffHBSM, ffHSM,
                                                      *sum(
                                                           ([(ggHg10gi2[ch], 1), (ttHg10gi2[ch], 1)]
                                                               for ch in channels),
                                                            []
                                                           )
                                                     )
            BSMVVH = self.ComponentTemplateSumInGroup("", VVHBSM, VVHSM,
                                                      *sum(
                                                           ([(VBFg10gi4[ch], 1), (ZHg10gi4[ch], 1), (WHg10gi4[ch], 1)]
                                                               for ch in channels),
                                                            []
                                                           )
                                                     )
            mix_pffH = self.ComponentTemplateSumInGroup("", ffHmix_p, ffHSM,
                                                        *sum(
                                                             ([(ggHg12gi0[ch], g1_mix**2), (ggHg10gi2[ch], (gi_mix/gi_ggHBSM)**2), (ggHg11gi1[ch],  g1_mix*gi_mix/gi_ggHBSM),
                                                               (ttHg12gi0[ch], g1_mix**2), (ttHg10gi2[ch], (gi_mix/gi_ggHBSM)**2), (ttHg11gi1[ch],  g1_mix*gi_mix/gi_ggHBSM)]
                                                                for ch in channels),
                                                              []
                                                             )
                                                       )
            mix_pVVH = self.ComponentTemplateSumInGroup("", VVHmix_p, VVHSM,
                                                        *sum(
                                                             ([(template, g1_mix**(4-j) * (+gi_mix)**j) for j, template in enumerate(VBFpieces[ch])]
                                                            + [(template, g1_mix**(4-j) * (+gi_mix)**j) for j, template in enumerate(ZHpieces[ch])]
                                                            + [(template, g1_mix**(4-j) * (+gi_mix)**j) for j, template in enumerate(WHpieces[ch])]
                                                                for ch in channels),
                                                              []
                                                             )
                                                       )
            mix_mffH = self.ComponentTemplateSumInGroup("", ffHmix_m, ffHSM,
                                                        *sum(
                                                             ([(ggHg12gi0[ch], g1_mix**2), (ggHg10gi2[ch], (gi_mix/gi_ggHBSM)**2), (ggHg11gi1[ch], -g1_mix*gi_mix/gi_ggHBSM),
                                                               (ttHg12gi0[ch], g1_mix**2), (ttHg10gi2[ch], (gi_mix/gi_ggHBSM)**2), (ttHg11gi1[ch], -g1_mix*gi_mix/gi_ggHBSM)]
                                                                for ch in channels),
                                                              []
                                                             )
                                                       )
            mix_mVVH = self.ComponentTemplateSumInGroup("", VVHmix_m, VVHSM,
                                                        *sum(
                                                             ([(template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(VBFpieces[ch])]
                                                            + [(template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(ZHpieces[ch])]
                                                            + [(template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(WHpieces[ch])]
                                                                for ch in channels),
                                                              []
                                                             )
                                                       )

            if category == "Untagged":
                SMbottom = SMVVH
                BSMbottom = BSMVVH
                mix_pbottom = mix_pVVH
                mix_mbottom = mix_mVVH
                bottomtitle = "VBF+VH"
                bottomcolor = 4
                bottommu = muV
                toptitle = "ggH+t#bar{t}H"
                topcolor = ROOT.kOrange+10
            elif category in ("VBFtagged", "VHHadrtagged"):
                SMbottom = SMffH
                BSMbottom = BSMffH
                mix_pbottom = mix_pffH
                mix_mbottom = mix_mffH
                bottomtitle = "ggH+t#bar{t}H"
                bottomcolor = ROOT.kOrange+10
                bottommu = muf
                toptitle = "VBF+VH"
                topcolor = 4

            SMbottom = self.TemplateSum("{} SM".format(bottomtitle),
                                       (SMbottom, bottommu), (ZZ, 1), #which already has ZX
                                       linecolor=bottomcolor, linewidth=2)
            BSMbottom = self.TemplateSum("{} {} = 1".format(bottomtitle, self.analysis.title()),
                                         (BSMbottom, bottommu), (ZZ, 1),
                                         linecolor=bottomcolor, linewidth=2, linestyle=2)
            mix_pbottom = self.TemplateSum("{} {} = #plus 0.5".format(bottomtitle, self.analysis.title(superscript=superscript)),
                                           (mix_pbottom, bottommu), (ZZ, 1),
                                           linecolor=bottomcolor, linewidth=2, linestyle=2)
            mix_mbottom = self.TemplateSum("{} {} = #minus 0.5".format(bottomtitle, self.analysis.title(superscript=superscript)),
                                           (mix_mbottom, bottommu), (ZZ, 1),
                                           linecolor=bottomcolor, linewidth=2, linestyle=2)

            SM = self.TemplateSum("{} SM".format(toptitle),
                                  (SMffH, muf), (SMVVH, muV), (ZZ, 1), #which already has ZX
                                  linecolor=topcolor, linewidth=2)
            BSM = self.TemplateSum("{} {} = 1".format(toptitle, self.analysis.title()),
                                   (BSMffH, muf), (BSMVVH, muV), (ZZ, 1),
                                   linecolor=topcolor, linewidth=2, linestyle=2)
            mix_p = self.TemplateSum("{} {} = #plus 0.5".format(toptitle, self.analysis.title(superscript=superscript)),
                                     (mix_pffH, muf), (mix_pVVH, muV), (ZZ, 1),
                                     linecolor=topcolor, linewidth=2, linestyle=2)
            mix_m = self.TemplateSum("{} {} = #minus 0.5".format(toptitle, self.analysis.title(superscript=superscript)),
                                     (mix_mffH, muf), (mix_mVVH, muV), (ZZ, 1),
                                     linecolor=topcolor, linewidth=2, linestyle=2)
            templates[0:0] = [SMbottom, SM, BSMbottom, BSM, mix_pbottom, mix_p, mix_mbottom, mix_m] #will remove some later, depending on the discriminant

            if category in ("VBFtagged", "VHHadrtagged"):
                rebin = 5
        else: #animation
            ffHintegral  = self.IntegralSum(*sum(
                                                 ([(allggHg12gi0[ca,ch], g1_custom**2), (allggHg10gi2[ca,ch], (gi_custom/gi_ggHBSM)**2), (allggHg11gi1[ca,ch],  g1_custom*gi_custom/gi_ggHBSM),
                                                   (allttHg12gi0[ca,ch], g1_custom**2), (allttHg10gi2[ca,ch], (gi_custom/gi_ggHBSM)**2), (allttHg11gi1[ca,ch],  g1_custom*gi_custom/gi_ggHBSM)]
                                                    for ca, ch in itertools.product(categories, channels)),
                                                  []
                                                 ))
            VVHintegral  = self.IntegralSum(*sum((
                                                   [(template, g1_custom**(4-j) * (+gi_custom)**j) for j, template in enumerate(allVBFpieces[ca,ch])]
                                                 + [(template, g1_custom**(4-j) * (+gi_custom)**j) for j, template in enumerate(allZHpieces[ca,ch])]
                                                 + [(template, g1_custom**(4-j) * (+gi_custom)**j) for j, template in enumerate(allWHpieces[ca,ch])]
                                                    for ca, ch in itertools.product(categories, channels)),
                                                 []
                                                ))
            ffH = self.ComponentTemplateSumInGroup("", ffHintegral, ffHSM,
                                                   *sum(
                                                        ([(ggHg12gi0[ch], g1_custom**2), (ggHg10gi2[ch], (gi_custom/gi_ggHBSM)**2), (ggHg11gi1[ch],  g1_custom*gi_custom/gi_ggHBSM),
                                                          (ttHg12gi0[ch], g1_custom**2), (ttHg10gi2[ch], (gi_custom/gi_ggHBSM)**2), (ttHg11gi1[ch],  g1_custom*gi_custom/gi_ggHBSM)]
                                                           for ch in channels),
                                                         []
                                                        )
                                                  )
            VVH = self.ComponentTemplateSumInGroup("", VVHintegral, VVHSM,
                                                   *sum(
                                                        ([(template, g1_custom**(4-j) * (+gi_custom)**j) for j, template in enumerate(VBFpieces[ch])]
                                                       + [(template, g1_custom**(4-j) * (+gi_custom)**j) for j, template in enumerate(ZHpieces[ch])]
                                                       + [(template, g1_custom**(4-j) * (+gi_custom)**j) for j, template in enumerate(WHpieces[ch])]
                                                           for ch in channels),
                                                         []
                                                        )
                                                  )
            if category == "Untagged":
                bottom = VVH
                bottomtitle = "VBF+VH"
                bottomcolor = 4
                bottommu = muV
                toptitle = "ggH+t#bar{t}H"
                topcolor = ROOT.kOrange+10
            elif category in ("VBFtagged", "VHHadrtagged"):
                bottom = ffH
                bottomtitle = "ggH+t#bar{t}H"
                bottomcolor = ROOT.kOrange+10
                bottommu = muf
                toptitle = "VBF+VH"
                topcolor = 4

            bottom = self.TemplateSum(bottomtitle,
                                     (bottom, bottommu), (ZZ, 1), #which already has ZX
                                     linecolor=bottomcolor, linewidth=2)
            top = self.TemplateSum(toptitle,
                                   (VVH, muV), (ffH, muf), (ZZ, 1), #which already has ZX
                                   linecolor=topcolor, linewidth=2)

            templates[0:0] = [bottom, top] #will remove some later, depending on the discriminant

            if category in ("VBFtagged", "VHHadrtagged"):
                rebin = 5

        if self.enrichstatus == "impoverish" and config.showblinddistributions or config.unblinddistributions:
            data = self.TemplateSum("",
                                    *((self.TemplateFromFile(
                                                             category, self.enrichstatus, self.normalization,
                                                             self.analysis, ch, "data", self.production, linecolor=1
                                                            ), 1) for ch in channels)
                                   )
            data = [style.asymmerrorsfromhistogram(data.Projection(i, rebin=rebin), showemptyerrors=False) for i in range(3)]
            for g in data:
                g.SetLineColor(1)
                g.SetMarkerColor(1)
                g.SetLineStyle(1)
                g.SetLineWidth(1)
                g.SetMarkerStyle(20)
                g.SetMarkerSize(1.2)
            otherthingstodraw.append((data, "P"))
            otherthingstoaddtolegend.append((data, "Observed", "ep"))


    if not justoneproductionmode and not animation and not nicestyle:
        if self.enrichstatus == "impoverish" and config.showblinddistributions or config.unblinddistributions:
            templates += [
                          self.TemplateFromFile(category, self.enrichstatus, self.normalization, self.analysis, channel, "data", self.production)
                         ]

    for i, discriminant in enumerate(self.discriminants(category)):
        usetemplates = templates
        if nicestyle and not animation:
            usetemplates = templates[:]
            tf = TemplatesFile("2e2mu", "ggh", category, self.production, self.analysis)
            if discriminant == tf.mixdiscriminant and self.analysis == "fa3" or discriminant in (tf.mixdiscriminant, tf.purediscriminant) and self.analysis == "fL1":
                usetemplates.remove(BSM)
                usetemplates.remove(mix_m)
                usetemplates.remove(BSMbottom)
                usetemplates.remove(mix_mbottom)
            elif discriminant in (tf.mixdiscriminant, tf.purediscriminant) and self.analysis == "fa2":
                usetemplates.remove(BSM)
                usetemplates.remove(mix_p)
                usetemplates.remove(BSMbottom)
                usetemplates.remove(mix_pbottom)
            else:
                usetemplates.remove(mix_p)
                usetemplates.remove(mix_m)
                usetemplates.remove(mix_pbottom)
                usetemplates.remove(mix_mbottom)
        hstack = ROOT.THStack("{}{}".format(discriminant.name, saveasappend), "")
        legend = ROOT.TLegend(*legendargs[i])
        style.applylegendstyle(legend)

        for thing in otherthingstoaddtolegend:
            if isinstance(thing[0], collections.Sequence) and len(thing[0]) == 3: thing = [thing[0][i]] + list(thing[1:])
            legend.AddEntry(*thing)

        for template in usetemplates:
            if template.linecolor:
                if floor:
                    template.Floor()
                if template.title:
                    template.AddToLegend(legend)
                hstack.Add(template.Projection(i, rebin=rebin), template.hstackoption)

        if nicestyle:
            ymax = style.ymax((hstack, "nostack"), (data[i], "P"))
            hstack.SetMaximum(ymax*1.25)

        c1.cd()
        hstack.Draw("nostack")
        hstack.GetXaxis().SetTitle(discriminant.title)
        hstack.GetYaxis().SetTitle(
                                   "Events / {:.2f}".format(
                                                            (hstack.GetXaxis().GetXmax() - hstack.GetXaxis().GetXmin()) / hstack.GetXaxis().GetNbins()
                                                           )
                                  )
        if nicestyle and discriminant.name == "D_CP_decay": hstack.GetXaxis().SetRangeUser(-.4, .4)
        style.applyaxesstyle(hstack)
        if nicestyle:
            style.cuttext(self.enrichstatus.cuttext(), x1=.48+.03*(discriminant.name=="D_bkg"), x2=.58+.03*(discriminant.name=="D_bkg"))
            style.CMS("Preliminary", float(Luminosity("fordata", self.production)))
        legend.Draw()
        for thing, option in otherthingstodraw:
            if isinstance(thing, collections.Sequence) and len(thing) == 3: thing = thing[i]
            if thing:
                thing.Draw(option)
        try:
            os.makedirs(os.path.join(saveasdir, subdir))
        except OSError:
            pass
        for ext in exts:
            c1.SaveAs(os.path.join(saveasdir, subdir, "{}{}.{}".format(discriminant.name, saveasappend, ext)))

  def saveasdir(self, *categoryandchannel):
      categoryandchannel = self.CategoryAndChannel(*categoryandchannel)
      assert self.normalization == "rescalemixtures"
      return os.path.join(config.plotsbasedir, "templateprojections", self.enrichstatus.dirname(), "{}_{}/{}/{}".format(self.analysis, self.production, categoryandchannel.category, categoryandchannel.channel))
  def saveasdir_niceplots(self, category):
      assert self.normalization == "rescalemixtures" and len(config.productionsforcombine) == 1
      return os.path.join(config.plotsbasedir, "templateprojections", "niceplots_new", self.enrichstatus.dirname(), "{}/{}".format(self.analysis, Category(category)))

  def discriminants(self, category):
      return TemplatesFile("2e2mu", self.shapesystematic, "ggh", self.analysis, self.production, category).discriminants

  class AnimationStep(object):
    def __init__(self, productionmodeforfai, analysis, fai, delay, muV=1, muf=1, deltaNLL=None):
        self.analysis = analysis
        self.delay = delay
        self.sample = samplewithfai("ggH", self.analysis, fai, productionmodeforfai=productionmodeforfai)
        self.fai_decay = self.sample.fai("ggH", self.analysis)
        self.muV = muV
        self.muf = muf
        self.deltaNLL = deltaNLL
    def __cmp__(self, other):
        assert self.analysis == other.analysis
        return cmp((self.fai_decay, self.delay), (other.fai_decay, other.delay))
    def fai(self, *args, **kwargs):
        return self.sample.fai(*args, **kwargs)
    @property
    @cache
    def excludedtext(self):
        allowed = run1info.isallowed2sigma(self.analysis, self.fai_decay)
        if allowed: return
        x1, y1, x2, y2 = .2, .7, .5, .8
        pt = ROOT.TPaveText(x1, y1, x2, y2, "brNDC")
        pt.SetBorderSize(0)
        pt.SetFillStyle(0)
        pt.SetTextAlign(12)
        pt.SetTextFont(42)
        pt.SetTextSize(0.045)
        text = pt.AddText("Excluded in Run 1")
        text.SetTextColor(2)
        return pt

    @property
    @cache
    def niceplotsinfo(self):
        x1, y1, x2, y2 = .65, .57, .85, .9
        pt = ROOT.TPaveText(x1, y1, x2, y2, "brNDC")
        pt.SetBorderSize(0)
        pt.SetFillStyle(0)
        pt.SetTextAlign(12)
        pt.SetTextFont(42)
        pt.SetTextSize(0.045)
        pt.AddText("{}={:.2f}".format(self.analysis.title(superscript="dec"), self.fai_decay))
        pt.AddText("{}={:.2f}".format(self.analysis.title(superscript="VBF"), self.fai("VBF", self.analysis)))
        pt.AddText("{}={:.2f}".format(self.analysis.title(superscript="VH"), self.fai("VH", self.analysis)))
        pt.AddText("{}={:.2f}".format("#mu_{V}", self.muV))
        pt.AddText("{}={:.2f}".format("#mu_{f}", self.muf))
        pt.AddText("{}={:.2f}".format("-2#Deltaln L", self.deltaNLL))
        return pt

  def animation(self, category, channel, floor=False, scantree=None):
      tmpdir = mkdtemp()

      nsteps = 200
      animation = []
      nicestyle = (scantree is not None)

      finaldir = os.path.join(self.saveasdir(category, channel), "animation")

      if nicestyle:
          if channel != "2e2mu": return

          finaldir = os.path.join(self.saveasdir_niceplots(category), "animation")

          minNLL, faiforminNLL = min((scantree.deltaNLL+scantree.nll+scantree.nll0, scantree.CMS_zz4l_fai1) for entry in scantree)
          VBFmix = samplewithfai("ggH", self.analysis, 0.5, productionmodeforfai="VBF").fai("ggH", self.analysis)
          VHmix = samplewithfai("ggH", self.analysis, 0.5, productionmodeforfai="VH").fai("ggH", self.analysis)
          pauses = [
                    min((scantree.CMS_zz4l_fai1 for entry in scantree), key = lambda x: abs(x-_))
                       for _ in (0, -1, 1, .5, -.5, VBFmix, -VBFmix, VHmix, -VHmix) if _ in (0, 1, -1)
                   ]
          for entry in scantree:
              animation.append(
                               self.AnimationStep(
                                                  "ggH",
                                                  self.analysis,
                                                  scantree.CMS_zz4l_fai1,
                                                  100 if scantree.CMS_zz4l_fai1 == faiforminNLL or scantree.CMS_zz4l_fai1 in pauses
                                                     else 10,
                                                  muV = scantree.r_VVH,
                                                  muf = scantree.r_ffH,
                                                  deltaNLL = 2*(scantree.deltaNLL+scantree.nll+scantree.nll0 - minNLL)
                                                 )
                              )
      else:
          for productionmode in "VBF", "VH", "ggH":
              animation += [
                            self.AnimationStep(
                                               productionmode,
                                               self.analysis,
                                               -1+(2.*i/nsteps),
                                               50 if ((-1+(2.*i/nsteps))*2).is_integer() else 5
                                              ) for i in range(nsteps+1)
                           ]

      animation = sorted(set(animation))
      while True:
          for step, nextstep in pairwise(animation[:]):
              if step.fai_decay == nextstep.fai_decay:
                  animation.remove(step)
                  break
          else:
              break

      kwargs_base = {
                     "saveasdir": tmpdir,
                     "floor": floor,
                     "nicestyle": nicestyle,
                    }

      for i, step in enumerate(animation):
          kwargs = kwargs_base.copy()
          kwargs["customfaiforanimation"] = step.fai_decay, "ggH"
          kwargs["saveasappend"] = i
          if nicestyle:
              kwargs["otherthingstodraw"] = [(step.niceplotsinfo, "")]
          else:
              kwargs["otherthingstodraw"] = [(step.excludedtext, "")]
          kwargs["muVmuf"] = step.muV, step.muf
          self.projections(category, channel, **kwargs)

      for discriminant in self.discriminants(category):
          convertcommand = ["gm", "convert", "-loop", "0"]
          lastdelay = None
          for i, step in enumerate(animation):
              if step.delay != lastdelay:
                  convertcommand += ["-delay", str(step.delay)]
              convertcommand.append(os.path.join(tmpdir, "{}{}.gif".format(discriminant.name, i)))
          try:
              os.makedirs(finaldir)
          except OSError:
              pass
          convertcommand.append(os.path.join(finaldir, "{}.gif".format(discriminant.name)))
          subprocess.check_call(convertcommand)

def projections(*args):
    Projections(*args).projections()

if __name__ == "__main__":
  def projections():
#    yield Projections("170203", "2e2mu", "fa3", "rescalemixtures", "fullrange", "VHHadrtagged")
#    return
    for production in config.productionsforcombine:
      for analysis in analyses:
        for normalization in normalizations:
          for enrichstatus in enrichstatuses:
            if normalization != "rescalemixtures": continue   #uncomment this to get the niceplots fast
            yield Projections(analysis, normalization, production, enrichstatus)

  length = len(list(projections()))
  for i, (p, ch, ca) in enumerate(itertools.product(projections(), channels, categories), start=1):
    scantree = ROOT.TChain("limit")
    for filename in glob.glob(os.path.join(config.repositorydir, "CMSSW_7_6_5", "src", "HiggsAnalysis", "HZZ4l_Combination",
                                       "CreateDatacards", "cards_{}_Feb28_mu".format(p.analysis), "higgsCombine_obs_*.root")):
        if re.match("higgsCombine_obs_lumi[0-9.]*(_[0-9]*,[-.0-9]*,[-.0-9]*)*.MultiDimFit.mH125.root", os.path.basename(filename)):
            scantree.Add(filename)
    #p.projections(ch, ca, nicestyle=True)
    p.animation(ca, ch, scantree=scantree)
    #p.projections(ch, ca)
    #p.projections(ch, ca, subdir="ggH", productionmode="ggH")
    #p.projections(ch, ca, subdir="VBF", productionmode="VBF")
    #p.projections(ch, ca, subdir="VH",  productionmode="VH")
    #if p.enrichstatus == "fullrange":
    #  p.animation(ca, ch)
    print i, "/", length*len(channels)*len(categories)
