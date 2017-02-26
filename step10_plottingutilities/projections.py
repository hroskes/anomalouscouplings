#!/usr/bin/env python
import abc
import collections
from helperstuff import config, constants, run1info
from helperstuff.combinehelpers import getrate, gettemplate
from helperstuff.enums import analyses, Analysis, categories, Category, Channel, channels, EnumItem, MultiEnum, MultiEnumABCMeta, MyEnum, Production, ProductionMode, productions, ShapeSystematic
from helperstuff.samples import ReweightingSample, samplewithfai
import helperstuff.rootoverloads.histogramaxisnumbers, helperstuff.rootoverloads.histogramfloor
import helperstuff.style
from helperstuff.templates import IntTemplate, Template, TemplatesFile
from helperstuff.utilities import cache, tfiles, pairwise
import itertools
from math import copysign, sqrt
import os
import ROOT
import subprocess
from tempfile import mkdtemp

class TemplateForProjection(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, *args, **kwargs):
        self.__initedhistogram = False
        self.projections = {}
        super(TemplateForProjection, self).__init__(*args, **kwargs)

    @abc.abstractmethod
    def inithistogram(self):
        self.__initedhistogram = True
        for attr in self.histogramattrs:
            assert hasattr(self, attr)

    def Projection(self, i):
        if i not in self.projections:
            result = self.h.Projection(i)
            result.SetName(self.discriminants[i].name)
            result.SetTitle(self.discriminants[i].title)
            result.SetXTitle(self.discriminants[i].title)
            result.SetLineColor(self.color)
            self.projections[i] = result
        return self.projections[i]

    def AddToLegend(self, legend):
        legend.AddEntry(self.Projection(0), self.title, "l")

    def Scale(self, *args, **kwargs):
        oldintegral = self.h.Integral()
        return self.h.Scale(*args, **kwargs)
        newintegral = self.h.Integral()
        self.integral *= newintegral / oldintegral

    def Integral(self):
        return self.integral

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
        return "h", "hstackoption", "integral"

class Normalization(MyEnum):
    enumname = "normalization"
    enumitems = [
                 EnumItem("rescalemixtures"),
                 EnumItem("areanormalize"),
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
    def __init__(self, color, *args, **kwargs):
        super(BaseTemplateFromFile, self).__init__(*args, **kwargs)
        self.color = color
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
            self.h.SetMarkerColor(self.color)
            self.h.SetMarkerStyle(20)
        else:
            self.hstackoption = "hist"

        self.Scale(scalefactor)

        self.integral = self.h.Integral()

        self.doenrich()

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
    def __init__(self, title, color, *templatesandfactors):
        self.title = title
        self.color = color
        self.templatesandfactors = templatesandfactors
        analyses = {t[0].analysis for t in templatesandfactors}
        assert len(analyses) == 1
        self.analysis = analyses.pop()
        super(TemplateSumBase, self).__init__()

    def __hash__(self):
        return hash((self.title, self.color, self.templatesandfactors))

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
        self.h = None
        self.integral = 0
        for template, factor in self.templatesandfactors:
            if not factor: continue
            self.integral += template.integral*factor
            if self.h is None:
                self.h = template.h.Clone()
                self.h.Scale(factor)
            else:
                self.h.Add(template.h, factor)
        self.hstackoption = "hist"
        super(TemplateSum, self).inithistogram()

class IntegralSum(TemplateSumBase):
    counter = 0
    def __init__(self, *templatesandfactors):
        super(IntegralSum, self).__init__("", 0, *templatesandfactors)
    @classmethod
    def hname(cls):
        cls.counter += 1
        return "h{}".format(cls.counter)
    def inithistogram(self):
        self.hstackoption = None
        self.integral = 0
        for template, factor in self.templatesandfactors:
            self.integral += template.Integral() * factor
        self.h = ROOT.TH1D(self.hname(), "", 1, 0, 1)
        self.h.Fill(.5, self.integral)
        assert abs(self.h.Integral() - self.integral)/self.integral < 1e-10, (self.h.Integral(), self.integral)
        super(IntegralSum, self).inithistogram()

class ComponentTemplateSum(TemplateSum):
    def __init__(self, title, color, SMintegral, *templatesandfactors):
        """
        Works the same as TemplateSum, but rescales
        """
        self.SMintegral = SMintegral
        super(ComponentTemplateSum, self).__init__(title, color, *templatesandfactors)
        for template, factor in templatesandfactors:
            if not isinstance(template, BaseTemplateFromFile):
                raise TypeError("{} can only come from TemplatesFromFiles".format(type(self).__name__))

    def __hash__(self):
        return hash((self.SMintegral, super(ComponentTemplateSum, self).__hash__()))

    def inithistogram(self):
        super(ComponentTemplateSum, self).inithistogram()
        if isinstance(self.SMintegral, TemplateForProjection): self.SMintegral = self.SMintegral.Integral()
        if self.Integral():
            self.Scale(self.SMintegral / self.Integral())

class ComponentTemplateSumInGroup(TemplateSum):
    def __init__(self, title, color, mygroup, SMgroup, *templatesandfactors):
        self.mygroup = mygroup
        self.SMgroup = SMgroup
        super(ComponentTemplateSumInGroup, self).__init__(title, color, *templatesandfactors)
        for template, factor in templatesandfactors:
            if not isinstance(template, BaseTemplateFromFile):
                raise TypeError("{} can only come from TemplatesFromFiles".format(type(self).__name__))

    def __hash__(self):
        return hash((self.mygroup, self.SMgroup, super(ComponentTemplateSumInGroup, self).__hash__()))

    def inithistogram(self):
        super(ComponentTemplateSumInGroup, self).inithistogram()
        if self.Integral():
            self.Scale(self.SMgroup.Integral() / self.mygroup.Integral())

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
  def TemplateFromFile(self, *args):
    return TemplateFromFile(*args)
  @cache
  def IntTemplateFromFile(self, *args):
    return IntTemplateFromFile(*args)
  @cache
  def TemplateSum(self, *args):
    return TemplateSum(*args)
  @cache
  def IntegralSum(self, *args):
    return IntegralSum(*args)
  @cache
  def ComponentTemplateSum(self, *args):
    return ComponentTemplateSum(*args)
  @cache
  def ComponentTemplateSumInGroup(self, *args):
    return ComponentTemplateSumInGroup(*args)

  def projections(self, *categoryandchannel, **kwargs):
    class _(MultiEnum): enums = (Category, Channel)
    info = _(*categoryandchannel)
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
    floor = False
    for kw, kwarg in kwargs.iteritems():
       if kw == "ggHfactor":
           ggHfactor = float(kwarg)
       elif kw == "VBFfactor":
           VBFfactor = float(kwarg)
       elif kw == "VHfactor":
           VHfactor = float(kwarg)
       elif kw == "ttHfactor":
           ttHfactor = float(kwarg)
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

    for ca, ch in itertools.product(categories, channels):
      ggHg12gi0[ca,ch] = self.TemplateFromFile(   0, "ggH", ca, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
      ggHg10gi2[ca,ch] = self.TemplateFromFile(   0, "ggH", ca, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, BSMhypothesis)
      ggHg11gi1[ca,ch] = self.IntTemplateFromFile(0, "ggH", ca, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g11gi1")

      ttHg12gi0[ca,ch] = self.TemplateFromFile(   0, "ttH", ca, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
      ttHg10gi2[ca,ch] = self.TemplateFromFile(   0, "ttH", ca, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, BSMhypothesis)
      ttHg11gi1[ca,ch] = self.IntTemplateFromFile(0, "ttH", ca, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g11gi1")

      VBFg14gi0[ca,ch] = self.TemplateFromFile(   0, "VBF", ca, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
      VBFg13gi1[ca,ch] = self.IntTemplateFromFile(0, "VBF", ca, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g13gi1")
      VBFg12gi2[ca,ch] = self.IntTemplateFromFile(0, "VBF", ca, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g12gi2")
      VBFg11gi3[ca,ch] = self.IntTemplateFromFile(0, "VBF", ca, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g11gi3")
      VBFg10gi4[ca,ch] = self.TemplateFromFile(   0, "VBF", ca, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, BSMhypothesis)

      VBFpieces[ca,ch] = [VBFg14gi0[ca,ch], VBFg13gi1[ca,ch], VBFg12gi2[ca,ch], VBFg11gi3[ca,ch], VBFg10gi4[ca,ch]]

      ZHg14gi0[ca,ch] = self.TemplateFromFile(   0, "ZH", ca, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
      ZHg13gi1[ca,ch] = self.IntTemplateFromFile(0, "ZH", ca, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g13gi1")
      ZHg12gi2[ca,ch] = self.IntTemplateFromFile(0, "ZH", ca, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g12gi2")
      ZHg11gi3[ca,ch] = self.IntTemplateFromFile(0, "ZH", ca, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g11gi3")
      ZHg10gi4[ca,ch] = self.TemplateFromFile(   0, "ZH", ca, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, BSMhypothesis)

      ZHpieces[ca,ch] = [ZHg14gi0[ca,ch], ZHg13gi1[ca,ch], ZHg12gi2[ca,ch], ZHg11gi3[ca,ch], ZHg10gi4[ca,ch]]

      WHg14gi0[ca,ch] = self.TemplateFromFile(   0, "WH", ca, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
      WHg13gi1[ca,ch] = self.IntTemplateFromFile(0, "WH", ca, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g13gi1")
      WHg12gi2[ca,ch] = self.IntTemplateFromFile(0, "WH", ca, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g12gi2")
      WHg11gi3[ca,ch] = self.IntTemplateFromFile(0, "WH", ca, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, "g11gi3")
      WHg10gi4[ca,ch] = self.TemplateFromFile(   0, "WH", ca, self.enrichstatus, self.normalization, self.production, ch, self.shapesystematic, self.analysis, BSMhypothesis)

      WHpieces[ca,ch] = [WHg14gi0[ca,ch], WHg13gi1[ca,ch], WHg12gi2[ca,ch], WHg11gi3[ca,ch], WHg10gi4[ca,ch]]

    ffHSM     = self.IntegralSum(*sum(
                                      ([(ggHg12gi0[ca,ch], 1), (ttHg12gi0[ca,ch], 1)]
                                          for ca,ch in itertools.product(categories, channels)),
                                       []
                                      ))
    ffHBSM    = self.IntegralSum(*sum(
                                      ([(ggHg10gi2[ca,ch], 1), (ttHg10gi2[ca,ch], 1)]
                                         for ca,ch in itertools.product(categories, channels)),
                                      []
                                     ))
    ffHmix_p  = self.IntegralSum(*sum(
                                      ([(ggHg12gi0[ca,ch], g1_mix**2), (ggHg10gi2[ca,ch], (gi_mix/gi_ggHBSM)**2), (ggHg11gi1[ca,ch],  g1_mix*gi_mix/gi_ggHBSM),
                                        (ttHg12gi0[ca,ch], g1_mix**2), (ttHg10gi2[ca,ch], (gi_mix/gi_ggHBSM)**2), (ttHg11gi1[ca,ch],  g1_mix*gi_mix/gi_ggHBSM)]
                                         for ca, ch in itertools.product(categories, channels)),
                                       []
                                      ))
    ffHmix_m  = self.IntegralSum(*sum(
                                             ([(ggHg12gi0[ca,ch], g1_mix**2), (ggHg10gi2[ca,ch], (gi_mix/gi_ggHBSM)**2), (ggHg11gi1[ca,ch], -g1_mix*gi_mix/gi_ggHBSM),
                                               (ttHg12gi0[ca,ch], g1_mix**2), (ttHg10gi2[ca,ch], (gi_mix/gi_ggHBSM)**2), (ttHg11gi1[ca,ch], -g1_mix*gi_mix/gi_ggHBSM)]
                                                for ca, ch in itertools.product(categories, channels)),
                                              []
                                            ))

    VVHSM     = self.IntegralSum(*sum(
                                             ([(VBFg14gi0[ca,ch], 1), (ZHg14gi0[ca,ch], 1), (WHg14gi0[ca,ch], 1)]
                                                for ca,ch in itertools.product(categories, channels)),
                                             []
                                            ))
    VVHBSM    = self.IntegralSum(*sum(
                                             ([(VBFg10gi4[ca,ch], 1), (ZHg10gi4[ca,ch], 1), (WHg10gi4[ca,ch], 1)]
                                                for ca,ch in itertools.product(categories, channels)),
                                             []
                                            ))
    VVHmix_p  = self.IntegralSum(*sum((
                                               [(template, g1_mix**(4-j) * (+gi_mix)**j) for j, template in enumerate(VBFpieces[ca,ch])]
                                             + [(template, g1_mix**(4-j) * (+gi_mix)**j) for j, template in enumerate(VBFpieces[ca,ch])]
                                             + [(template, g1_mix**(4-j) * (+gi_mix)**j) for j, template in enumerate(VBFpieces[ca,ch])]
                                                for ca,ch in itertools.product(categories, channels)),
                                             []
                                            ))
    VVHmix_m  = self.IntegralSum(*sum((
                                               [(template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(VBFpieces[ca,ch])]
                                             + [(template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(VBFpieces[ca,ch])]
                                             + [(template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(VBFpieces[ca,ch])]
                                                for ca,ch in itertools.product(categories, channels)),
                                             []
                                            ))

    ggHSM = self.TemplateSum("ggH SM", 0, (ggHg12gi0[category,channel], 1))
    ggHBSM = self.ComponentTemplateSumInGroup("ggH {}=1".format(fainame), 0, ffHBSM, ffHSM, (ggHg10gi2[category,channel], 1))
    ggHmix_p  = self.ComponentTemplateSumInGroup("ggH {}=#plus0.5" .format(fainame), 0, ffHmix_p, ffHSM,
                                                 (ggHg12gi0[category,channel], g1_mix**2), (ggHg10gi2[category,channel], (gi_mix/gi_ggHBSM)**2), (ggHg11gi1[category,channel],  g1_mix*gi_mix/gi_ggHBSM))
    ggHmix_m  = self.ComponentTemplateSumInGroup("ggH {}=#minus0.5".format(fainame), 0, ffHmix_m, ffHSM,
                                                 (ggHg12gi0[category,channel], g1_mix**2), (ggHg10gi2[category,channel], (gi_mix/gi_ggHBSM)**2), (ggHg11gi1[category,channel], -g1_mix*gi_mix/gi_ggHBSM))
    ttHSM = self.TemplateSum("ttH SM", 0, (ttHg12gi0[category,channel], 1))
    ttHBSM = self.ComponentTemplateSumInGroup("ttH {}=1".format(fainame), 0, ffHBSM, ffHSM, (ttHg10gi2[category,channel], 1))
    ttHmix_p  = self.ComponentTemplateSumInGroup("ttH {}=#plus0.5" .format(fainame), 0, ffHmix_p, ffHSM,
                                                 (ttHg12gi0[category,channel], g1_mix**2), (ttHg10gi2[category,channel], (gi_mix/gi_ggHBSM)**2), (ttHg11gi1[category,channel],  g1_mix*gi_mix/gi_ggHBSM))
    ttHmix_m  = self.ComponentTemplateSumInGroup("ttH {}=#minus0.5".format(fainame), 0, ffHmix_m, ffHSM,
                                                 (ttHg12gi0[category,channel], g1_mix**2), (ttHg10gi2[category,channel], (gi_mix/gi_ggHBSM)**2), (ttHg11gi1[category,channel], -g1_mix*gi_mix/gi_ggHBSM))

    VBFSM = self.TemplateSum("VBF SM", 0, (VBFg14gi0[category,channel], 1))
    VBFBSM = self.ComponentTemplateSumInGroup("VBF {}=1".format(fainame), 0, VVHBSM, VVHSM, (VBFg10gi4[category,channel], 1))
    VBFmix_p  = self.ComponentTemplateSumInGroup("VBF {}=#plus0.5" .format(fainame), 0, VVHmix_p, VVHSM,
                                            *((template, g1_mix**(4-j) * gi_mix**j) for j, template in enumerate(VBFpieces[category,channel]))
                                           )
    VBFmix_m  = self.ComponentTemplateSumInGroup("VBF {}=#minus0.5".format(fainame), 0, VVHmix_m, VVHSM,
                                            *((template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(VBFpieces[category, channel]))
                                           )

    VHSM = self.TemplateSum("VH SM", 0, (ZHg14gi0[category,channel], 1), (WHg14gi0[category,channel], 1))
    VHBSM = self.ComponentTemplateSumInGroup("VH {}=1".format(fainame), 0, VVHBSM, VVHSM, (ZHg10gi4[category,channel], 1), (WHg10gi4[category,channel], 1))
    VHmix_p  = self.ComponentTemplateSumInGroup("VH {}=#plus0.5" .format(fainame), 0, VVHmix_p, VVHSM,
                                            *(
                                                [(template, g1_mix**(4-j) * gi_mix**j) for j, template in enumerate(ZHpieces[category,channel])]
                                              + [(template, g1_mix**(4-j) * gi_mix**j) for j, template in enumerate(WHpieces[category,channel])]
                                             )
                                           )
    VHmix_m  = self.ComponentTemplateSumInGroup("VH {}=#minus0.5" .format(fainame), 0, VVHmix_m, VVHSM,
                                            *(
                                                [(template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(ZHpieces[category,channel])]
                                              + [(template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(WHpieces[category,channel])]
                                             )
                                           )

    SM    = self.TemplateSum("SM",                                 1,             (ggHSM,    ggHfactor), (VBFSM,    VBFfactor), (VHSM,    VHfactor), (ttHSM,    ttHfactor))
    BSM   = self.TemplateSum("{}=1".format(self.analysis.title()), ROOT.kCyan,    (ggHBSM,   ggHfactor), (VBFBSM,   VBFfactor), (VHBSM,   VHfactor), (ttHBSM,   ttHfactor))
    mix_p = self.TemplateSum("{}=#plus0.5".format(fainame),        ROOT.kGreen+3, (ggHmix_p, ggHfactor), (VBFmix_p, VBFfactor), (VHmix_p, VHfactor), (ttHmix_p, ttHfactor))
    mix_m = self.TemplateSum("{}=#minus0.5".format(fainame),       4,             (ggHmix_m, ggHfactor), (VBFmix_m, VBFfactor), (VHmix_m, VHfactor), (ttHmix_m, ttHfactor))

    qqZZ      = self.TemplateFromFile(6,              category, self.enrichstatus, self.normalization, self.analysis, channel, "qqZZ",    self.shapesystematic, self.production)
    ggZZ      = self.TemplateFromFile(ROOT.kOrange+6, category, self.enrichstatus, self.normalization, self.analysis, channel, "ggZZ",    self.shapesystematic, self.production)
    VBFbkg    = self.TemplateFromFile(ROOT.kViolet+7, category, self.enrichstatus, self.normalization, self.analysis, channel, "VBF bkg", self.shapesystematic, self.production)
    ZX        = self.TemplateFromFile(2,              category, self.enrichstatus, self.normalization, self.analysis, channel, "ZX",      self.shapesystematic, self.production)

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
        ggHcustom = self.ComponentTemplateSum("ggH ({}={}{:.2f})".format(self.analysis.title(superscript="dec"), plusminus, abs(customsample.fai("ggH", self.analysis))), 1, ggHSM.Integral(), (ggHSM, g1_custom**2), (ggHBSM, (gi_custom/gi_ggHBSM)**2), (ggHg11gi1,  g1_custom*gi_custom/gi_ggHBSM))
        VBFcustom = self.ComponentTemplateSum("VBF ({}={}{:.2f})".format(self.analysis.title(superscript="VBF"), plusminus, abs(customsample.fai("VBF", self.analysis))), 2, VBFSM.Integral(),
                                         *((template, g1_custom**(4-j) * gi_custom**j) for j, template in enumerate(VBFpieces))
                                        )
        VHcustom  = self.ComponentTemplateSum("VH ({}={}{:.2f})".format(self.analysis.title(superscript="VH"), plusminus, abs(customsample.fai("VH", self.analysis))), 4, ZHSM.Integral()+WHSM.Integral(),
                                         *(
                                             [(template, g1_custom**(4-j) * gi_custom**j) for j, template in enumerate(ZHpieces)]
                                           + [(template, g1_custom**(4-j) * gi_custom**j) for j, template in enumerate(WHpieces)]
                                          )
                                        )
        for t in ggHcustom, VBFcustom, VHcustom:
            t.Scale(1.0/t.Integral())

        templates += [
                      ggHcustom, VBFcustom, VHcustom
                     ]

    if ggHfactor == VBFfactor == VHfactor == 1 and not animation:
        templates += [
                      qqZZ, ggZZ, VBFbkg, ZX,
                     ]

    if self.enrichstatus == "impoverish" and config.showblinddistributions or config.unblinddistributions:
        templates += [
                      self.TemplateFromFile(1, category, self.enrichstatus, self.normalization, self.analysis, channel, "data", self.production)
                     ]

    c1 = ROOT.TCanvas()
    legend = ROOT.TLegend(.65, .6, .9, .9)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    for template in templates:
        if template.color:
            if self.normalization == "areanormalize":
                try:
                    template.Scale(1/template.Integral())
                except ZeroDivisionError:
                    pass
            if floor:
                template.Floor()
            if template.title:
                template.AddToLegend(legend)

    for i, discriminant in enumerate(self.discriminants(category)):
        hstack = ROOT.THStack("{}{}".format(discriminant.name, saveasappend), discriminant.title)
        for template in templates:
            if template.color:
                hstack.Add(template.Projection(i), template.hstackoption)
        hstack.Draw("nostack")
        hstack.GetXaxis().SetTitle(discriminant.title)
        legend.Draw()
        for thing in otherthingstodraw:
            if thing:
                thing.Draw()
        try:
            os.makedirs(os.path.join(saveasdir, subdir))
        except OSError:
            pass
        for ext in exts:
            c1.SaveAs(os.path.join(saveasdir, subdir, "{}{}.{}".format(discriminant.name, saveasappend, ext)))

  def saveasdir(self, categoryandchannel):
      assert self.normalization == "rescalemixtures"
      return os.path.join(config.plotsbasedir, "templateprojections", self.enrichstatus.dirname(), "{}_{}/{}/{}".format(self.analysis, self.production, categoryandchannel.category, categoryandchannel.channel))

  def discriminants(self, category):
      return TemplatesFile("2e2mu", self.shapesystematic, "ggh", self.analysis, self.production, category).discriminants

  class AnimationStep(object):
    def __init__(self, productionmodeforfai, analysis, fai, delay):
        self.analysis = analysis
        self.delay = delay
        sample = samplewithfai("ggH", self.analysis, fai, productionmodeforfai=productionmodeforfai)
        self.fai_decay = sample.fai("ggH", self.analysis)
    def __cmp__(self, other):
        assert self.analysis == other.analysis
        return cmp((self.fai_decay, self.delay), (other.fai_decay, other.delay))
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

  def animation(self, category, channel, floor=False):
      tmpdir = mkdtemp()

      nsteps = 200
      animation = []
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
                    }

      for i, step in enumerate(animation):
          kwargs = kwargs_base.copy()
          kwargs["customfaiforanimation"] = step.fai_decay, "ggH"
          kwargs["saveasappend"] = i
          kwargs["otherthingstodraw"] = [step.excludedtext]
          self.projections(category, channel, **kwargs)

      for discriminant in self.discriminants(category):
          convertcommand = ["gm", "convert", "-loop", "0"]
          lastdelay = None
          for i, step in enumerate(animation):
              if step.delay != lastdelay:
                  convertcommand += ["-delay", str(step.delay)]
              convertcommand.append(os.path.join(tmpdir, "{}{}.gif".format(discriminant.name, i)))
          try:
              os.makedirs(os.path.join(self.saveasdir, "animation"))
          except OSError:
              pass
          convertcommand.append(os.path.join(self.saveasdir, "animation", "{}.gif".format(discriminant.name)))
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
    p.projections(ch, ca)
    #p.projections(ch, ca, subdir="ggH", ggHfactor=1, VBFfactor=0, VHfactor=0)
    #p.projections(ch, ca, subdir="VBF", ggHfactor=0, VBFfactor=1, VHfactor=0)
    #p.projections(ch, ca, subdir="VH",  ggHfactor=0, VBFfactor=0, VHfactor=1)
    #if p.enrichstatus == "fullrange":
      #p.animation(ch, ca)
    print i, "/", length*len(channels)*len(categories)
