#!/usr/bin/env python
import abc
import collections
import glob
import itertools
from math import copysign, sqrt
import os
import pipes
import re
import shutil
import subprocess
import sys

import ROOT

from helperstuff import config, run1info, stylefunctions as style
from helperstuff.combinehelpers import getdatatree2015, getrate, gettemplate, Luminosity
from helperstuff.enums import analyses, Analysis, categories, Category, Channel, channels, EnumItem, MultiEnum, MultiEnumABCMeta, MyEnum, Production, ProductionMode, productions, ShapeSystematic
from helperstuff.samples import ReweightingSample, samplewithfai
from helperstuff.submitjob import submitjob
import helperstuff.rootoverloads.histogramaxisnumbers, helperstuff.rootoverloads.histogramfloor
from helperstuff.templates import IntTemplate, Template, TemplatesFile
from helperstuff.utilities import cache, tfiles, mkdtemp, pairwise

c1 = ROOT.TCanvas("cprojections", "", 8, 30, 800, 800)
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
        self.SetHistStyle(result)
        return result

    def SetHistStyle(self, h):
        h.SetLineColor(self.linecolor)
        h.SetLineStyle(self.linestyle)
        h.SetLineWidth(self.linewidth)
        h.SetFillColor(self.fillcolor)
        h.SetFillStyle(self.fillstyle)

    @cache
    def Project3D(self, option, rebinx=None, rebiny=None):
        if len(option) != 2 or not (set(option) <= set("xyz")):
            raise ValueError("Project3D only works for 2D projections")
        i = ["xyz".index(_) for _ in option]
        if rebinx is not None or rebiny is not None:
            result = self.Project3D(option).Clone()
            if rebinx is not None:
                result.RebinX(rebinx)
            if rebiny is not None:
                result.RebinY(rebiny)
            return result

        result = self.h.Project3D(option)
        result.SetName("_".join(self.discriminants[_].name for _ in i))
        result.SetTitle(" ".join(self.discriminants[_].title for _ in i))
        result.SetXTitle(self.discriminants[i[1]].title)
        result.SetYTitle(self.discriminants[i[0]].title)
        return result

    def AddToLegend(self, legend):
        legend.AddEntry(self.Projection(2), self.title, self.legendoption)

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
        return "h", "hstackoption", "cantakeintegral", "isDbkgonly"

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
        self.kwargs = kwargs
        for kw, kwarg in kwargs.items():
          if kw in ():
            setattr(self, kw, kwarg)
            del kwargs[kw]
        super(BaseTemplateFromFile, self).__init__(*args, **kwargs)
        self.didenrich = False
        #I don't know what the following line accomplishes
        # but it fixes the WH interference templates.
        #Otherwise they are wrong in some way that
        # I don't understand
        self.template.gettemplate()

    def __hash__(self):
        return hash((super(BaseTemplateFromFile, self).__hash__(), tuple(sorted(self.kwargs.iteritems()))))

    def inithistogram(self):
        self.isDbkgonly = False
        self.h = self.template.gettemplate().Clone()

        if self.productionmode in ["ggH", "VBF", "ZH", "WH"]:
            scalefactor = getrate(self.channel, self.category, self.productionmode, "fordata", self.production, self.analysis) / Template(self.production, self.category, self.analysis, self.channel, self.productionmode, self.analysis.purehypotheses[0]).gettemplate().Integral()
            numerator = getrate(self.channel, self.category, self.productionmode, "fordata", self.production, self.analysis)
            denominator = Template(self.production, self.category, self.analysis, self.channel, self.productionmode, self.analysis.purehypotheses[0]).gettemplate().Integral()
            scalefactor = numerator / denominator
        elif self.productionmode == "data" or self.h.Integral() == 0:
            scalefactor = 1
        elif self.h.Integral() == 0:
            scalefactor = 1
        else:
            numerator = getrate(self.channel, self.category, self.productionmode, "fordata", self.production, self.analysis)
            denominator = self.h.Integral()
            scalefactor = numerator / denominator

        if self.productionmode == "data":
            self.hstackoption = "P"
            self.h.SetMarkerColor(self.linecolor)
            self.h.SetMarkerStyle(20)
        else:
            self.hstackoption = "hist"

        self.h.Scale(scalefactor)

        self.doenrich()

        self.cantakeintegral = True#(self.enrichstatus == "fullrange")

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
        if not isinstance(self.title, basestring):
            raise TypeError("Title should be a string, not {}!!".format(self.title))
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

    @cache
    def Projection(self, i, rebin=None):
        if self.isDbkgonly and rebin is None:  #rebin is handled in super
            if i == 2:
                return self.h.Clone()
            else:
                raise ValueError("Can only get the 2nd (Z) projection of {}, not Projection({})".format(self, i))
        else:
            return super(TemplateSumBase, self).Projection(i, rebin=rebin)

class TemplateSum(TemplateSumBase):
    def inithistogram(self):
        self.cantakeintegral = True
        self.h = None

        if all(template.isDbkgonly for template, factor in self.templatesandfactors):
            self.isDbkgonly = True
        elif not any(template.isDbkgonly for template, factor in self.templatesandfactors):
            self.isDbkgonly = False
        else:
            raise ValueError("Some of the templates for {} are Dbkg only and some are not!".format(self))

        for template, factor in self.templatesandfactors:
            if not factor: continue
            if self.h is None:
                self.h = template.h.Clone()
                self.h.Scale(factor)
            else:
                self.h.Add(template.h, factor)
            self.cantakeintegral = self.cantakeintegral and template.cantakeintegral
        if self.isDbkgonly:
            self.SetHistStyle(self.h)
        self.hstackoption = "hist"
        super(TemplateSum, self).inithistogram()

class DbkgSum(TemplateSumBase):
    def inithistogram(self):
        self.cantakeintegral = True
        self.h = None
        self.isDbkgonly = True
        integral = 0
        for template, factor in self.templatesandfactors:
            if not factor: continue
            self.cantakeintegral = self.cantakeintegral and template.cantakeintegral
            if self.cantakeintegral:
                integral += template.Integral() * factor
            if self.h is None:
                self.h = template.Projection(2).Clone()
                self.h.Scale(factor)
            else:
                self.h.Add(template.Projection(2), factor)
        if self.cantakeintegral:
            assert abs(self.h.Integral() - integral)/integral < 1e-10, (self.h.Integral(), integral)
        self.SetHistStyle(self.h)
        self.hstackoption = "hist"
        super(DbkgSum, self).inithistogram()

class ComponentTemplateSumBase(TemplateSumBase):
    def __init__(self, title, SMintegral, *templatesandfactors, **kwargs):
        """
        Works the same as TemplateSum, but rescales
        """
        self.SMintegral = SMintegral
        super(ComponentTemplateSumBase, self).__init__(title, *templatesandfactors, **kwargs)
        for template, factor in templatesandfactors:
            if not isinstance(template, BaseTemplateFromFile):
                raise TypeError("{} can only come from TemplatesFromFiles".format(type(self).__name__))

    def __hash__(self):
        return hash((self.SMintegral, super(ComponentTemplateSumBase, self).__hash__()))

    def inithistogram(self):
        super(ComponentTemplateSumBase, self).inithistogram()
        if isinstance(self.SMintegral, TemplateForProjection):
            #if not self.cantakeintegral:
            #    raise ValueError("Have to be able to take integral of ComponentTemplateSum, or normalize to a number")
            self.SMintegral = self.SMintegral.Integral()
        if self.h.Integral():
            self.h.Scale(self.SMintegral / self.h.Integral())

class ComponentTemplateSum(ComponentTemplateSumBase, TemplateSum):  #the order is important!
    pass

class ComponentDbkgSum(ComponentTemplateSumBase, DbkgSum):  #the order is important!
    pass

class ComponentTemplateSumInGroupBase(TemplateSumBase):
    def __init__(self, title, mygroup, SMgroup, *templatesandfactors, **kwargs):
        self.mygroup = mygroup
        self.SMgroup = SMgroup
        super(ComponentTemplateSumInGroupBase, self).__init__(title, *templatesandfactors, **kwargs)
        for template, factor in templatesandfactors:
            if not isinstance(template, BaseTemplateFromFile):
                raise TypeError("{} can only come from TemplatesFromFiles".format(type(self).__name__))

    def __hash__(self):
        return hash((self.mygroup, self.SMgroup, super(ComponentTemplateSumInGroupBase, self).__hash__()))

    def inithistogram(self):
        super(ComponentTemplateSumInGroupBase, self).inithistogram()
        self.h.Scale(self.SMgroup.Integral() / self.mygroup.Integral())

class ComponentTemplateSumInGroup(ComponentTemplateSumInGroupBase, TemplateSum):  #the order is important!
    pass

class ComponentDbkgSumInGroup(ComponentTemplateSumInGroupBase, DbkgSum):  #the order is important!
    pass

class Projections(MultiEnum):
  enums = [Analysis, Normalization, ShapeSystematic, EnrichStatus]
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

  @classmethod
  @cache
  def TemplateFromFile(cls, *args, **kwargs):
    return TemplateFromFile(*args, **kwargs)
  @classmethod
  @cache
  def IntTemplateFromFile(cls, *args, **kwargs):
    return IntTemplateFromFile(*args, **kwargs)
  @classmethod
  @cache
  def TemplateSum(cls, *args, **kwargs):
    return TemplateSum(*args, **kwargs)
  @classmethod
  @cache
  def DbkgSum(cls, *args, **kwargs):
    return DbkgSum(*args, **kwargs)
  @classmethod
  @cache
  def ComponentTemplateSum(cls, *args, **kwargs):
    return ComponentTemplateSum(*args, **kwargs)
  @classmethod
  @cache
  def ComponentDbkgSum(cls, *args, **kwargs):
    return ComponentDbkgSum(*args, **kwargs)
  @classmethod
  @cache
  def ComponentTemplateSumInGroup(cls, *args, **kwargs):
    return ComponentTemplateSumInGroup(*args, **kwargs)
  @classmethod
  @cache
  def ComponentDbkgSumInGroup(cls, *args, **kwargs):
    return ComponentDbkgSumInGroup(*args, **kwargs)

  class CategoryAndChannel(MultiEnum): enums = (Category, Channel)

  def projections(self, *categoryandchannel, **kwargs):
    info = self.CategoryAndChannel(*categoryandchannel)
    category, channel = info.category, info.channel

    giname = self.analysis.couplingname
    BSMhypothesis = self.analysis.purehypotheses[1]

    ggHfactor = VBFfactor = VHfactor = ttHfactor = 1
    subdir = ""
    customfai = None
    animation = False
    saveasdir = None
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
    Dbkg_allcategories = False
    forWIN = False
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
           if nicestyle and channel != "2e2mu": return
       elif kw == "Dbkg_allcategories":
           Dbkg_allcategories = bool(int(kwarg))
           if Dbkg_allcategories and category != "Untagged": return
       elif kw == "forWIN":
           forWIN = kwarg
       else:
           raise TypeError("Unknown kwarg {}={}".format(kw, kwarg))

    if saveasdir is None:
        if Dbkg_allcategories:
            saveasdir = self.saveasdir_Dbkgsum(forWIN=forWIN)
        elif nicestyle:
            saveasdir = self.saveasdir_niceplots(category, forWIN=forWIN)
        else:
            saveasdir = self.saveasdir(info, forWIN=forWIN)

    if Dbkg_allcategories and not nicestyle: raise ValueError("Dbkg_allcategories requires nicestyle!")

    SMhypothesis = self.analysis.purehypotheses[0]
    gi_ggHBSM = getattr(ReweightingSample("ggH", BSMhypothesis), BSMhypothesis.couplingname)
    gi_VBFBSM = copysign((ReweightingSample("VBF", SMhypothesis).nominalJHUxsec / ReweightingSample("VBF", BSMhypothesis).nominalJHUxsec)**.25, gi_ggHBSM)
    gi_VHBSM = copysign(((ReweightingSample("WH", SMhypothesis).nominalJHUxsec + ReweightingSample("ZH", SMhypothesis).nominalJHUxsec) / (ReweightingSample("WH", BSMhypothesis).nominalJHUxsec + ReweightingSample("ZH", BSMhypothesis).nominalJHUxsec))**.25, gi_ggHBSM)
    gi_VVHBSM = copysign(((ReweightingSample("VBF", SMhypothesis).nominalJHUxsec + ReweightingSample("WH", SMhypothesis).nominalJHUxsec + ReweightingSample("ZH", SMhypothesis).nominalJHUxsec) / (ReweightingSample("VBF", BSMhypothesis).nominalJHUxsec + ReweightingSample("WH", BSMhypothesis).nominalJHUxsec + ReweightingSample("ZH", BSMhypothesis).nominalJHUxsec))**.25, gi_ggHBSM)
    if category == "UntaggedAC19":
        g1_mix = 1/sqrt(2)
        gi_mix = 1/sqrt(2)*gi_ggHBSM
        fainame = self.analysis.title(superscript="dec")
    elif category == "VBF2jTaggedAC19":
        g1_mix = 1/2**.25
        gi_mix = 1/2**.25 * gi_VBFBSM
        fainame = self.analysis.title(superscript="VBFdec")
    elif category == "VHHadrTaggedAC19":
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

    fullrange = EnrichStatus("fullrange")

    for ca, ch, production in itertools.product(categories, channels, config.productionsforcombine):
      allggHg12gi0[ca,ch,production] = self.TemplateFromFile(   "ggH", ca, fullrange, self.normalization, production, ch, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
      allggHg10gi2[ca,ch,production] = self.TemplateFromFile(   "ggH", ca, fullrange, self.normalization, production, ch, self.shapesystematic, self.analysis, BSMhypothesis)
      allggHg11gi1[ca,ch,production] = self.IntTemplateFromFile("ggH", ca, fullrange, self.normalization, production, ch, self.shapesystematic, self.analysis, "g11gi1")

      allttHg12gi0[ca,ch,production] = self.TemplateFromFile(   "ttH", ca, fullrange, self.normalization, production, ch, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
      allttHg10gi2[ca,ch,production] = self.TemplateFromFile(   "ttH", ca, fullrange, self.normalization, production, ch, self.shapesystematic, self.analysis, BSMhypothesis)
      allttHg11gi1[ca,ch,production] = self.IntTemplateFromFile("ttH", ca, fullrange, self.normalization, production, ch, self.shapesystematic, self.analysis, "g11gi1")

      allVBFg14gi0[ca,ch,production] = self.TemplateFromFile(   "VBF", ca, fullrange, self.normalization, production, ch, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
      allVBFg13gi1[ca,ch,production] = self.IntTemplateFromFile("VBF", ca, fullrange, self.normalization, production, ch, self.shapesystematic, self.analysis, "g13gi1")
      allVBFg12gi2[ca,ch,production] = self.IntTemplateFromFile("VBF", ca, fullrange, self.normalization, production, ch, self.shapesystematic, self.analysis, "g12gi2")
      allVBFg11gi3[ca,ch,production] = self.IntTemplateFromFile("VBF", ca, fullrange, self.normalization, production, ch, self.shapesystematic, self.analysis, "g11gi3")
      allVBFg10gi4[ca,ch,production] = self.TemplateFromFile(   "VBF", ca, fullrange, self.normalization, production, ch, self.shapesystematic, self.analysis, BSMhypothesis)

      allVBFpieces[ca,ch,production] = [allVBFg14gi0[ca,ch,production], allVBFg13gi1[ca,ch,production], allVBFg12gi2[ca,ch,production], allVBFg11gi3[ca,ch,production], allVBFg10gi4[ca,ch,production]]

      allZHg14gi0[ca,ch,production] = self.TemplateFromFile(   "ZH", ca, fullrange, self.normalization, production, ch, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
      allZHg13gi1[ca,ch,production] = self.IntTemplateFromFile("ZH", ca, fullrange, self.normalization, production, ch, self.shapesystematic, self.analysis, "g13gi1")
      allZHg12gi2[ca,ch,production] = self.IntTemplateFromFile("ZH", ca, fullrange, self.normalization, production, ch, self.shapesystematic, self.analysis, "g12gi2")
      allZHg11gi3[ca,ch,production] = self.IntTemplateFromFile("ZH", ca, fullrange, self.normalization, production, ch, self.shapesystematic, self.analysis, "g11gi3")
      allZHg10gi4[ca,ch,production] = self.TemplateFromFile(   "ZH", ca, fullrange, self.normalization, production, ch, self.shapesystematic, self.analysis, BSMhypothesis)

      allZHpieces[ca,ch,production] = [allZHg14gi0[ca,ch,production], allZHg13gi1[ca,ch,production], allZHg12gi2[ca,ch,production], allZHg11gi3[ca,ch,production], allZHg10gi4[ca,ch,production]]

      allWHg14gi0[ca,ch,production] = self.TemplateFromFile(   "WH", ca, fullrange, self.normalization, production, ch, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
      allWHg13gi1[ca,ch,production] = self.IntTemplateFromFile("WH", ca, fullrange, self.normalization, production, ch, self.shapesystematic, self.analysis, "g13gi1")
      allWHg12gi2[ca,ch,production] = self.IntTemplateFromFile("WH", ca, fullrange, self.normalization, production, ch, self.shapesystematic, self.analysis, "g12gi2")
      allWHg11gi3[ca,ch,production] = self.IntTemplateFromFile("WH", ca, fullrange, self.normalization, production, ch, self.shapesystematic, self.analysis, "g11gi3")
      allWHg10gi4[ca,ch,production] = self.TemplateFromFile(   "WH", ca, fullrange, self.normalization, production, ch, self.shapesystematic, self.analysis, BSMhypothesis)

      allWHpieces[ca,ch,production] = [allWHg14gi0[ca,ch,production], allWHg13gi1[ca,ch,production], allWHg12gi2[ca,ch,production], allWHg11gi3[ca,ch,production], allWHg10gi4[ca,ch,production]]

    ffHSM     = self.DbkgSum("ffHSM", *sum(
                                           ([(allggHg12gi0[ca,ch,production], 1), (allttHg12gi0[ca,ch,production], 1)]
                                              for ca, ch, production in itertools.product(categories, channels, config.productionsforcombine)),
                                            []
                                           ))
    ffHBSM    = self.DbkgSum("ffHBSM", *sum(
                                            ([(allggHg10gi2[ca,ch,production], 1), (allttHg10gi2[ca,ch,production], 1)]
                                               for ca, ch, production in itertools.product(categories, channels, config.productionsforcombine)),
                                             []
                                           ))
    ffHmix_p  = self.DbkgSum("ffHmix_p", *sum(
                                              ([(allggHg12gi0[ca,ch,production], g1_mix**2), (allggHg10gi2[ca,ch,production], (gi_mix/gi_ggHBSM)**2), (allggHg11gi1[ca,ch,production],  g1_mix*gi_mix/gi_ggHBSM),
                                                (allttHg12gi0[ca,ch,production], g1_mix**2), (allttHg10gi2[ca,ch,production], (gi_mix/gi_ggHBSM)**2), (allttHg11gi1[ca,ch,production],  g1_mix*gi_mix/gi_ggHBSM)]
                                                 for ca, ch, production in itertools.product(categories, channels, config.productionsforcombine)),
                                               []
                                              ))
    ffHmix_m  = self.DbkgSum("ffHmix_m", *sum(
                                              ([(allggHg12gi0[ca,ch,production], g1_mix**2), (allggHg10gi2[ca,ch,production], (gi_mix/gi_ggHBSM)**2), (allggHg11gi1[ca,ch,production], -g1_mix*gi_mix/gi_ggHBSM),
                                                (allttHg12gi0[ca,ch,production], g1_mix**2), (allttHg10gi2[ca,ch,production], (gi_mix/gi_ggHBSM)**2), (allttHg11gi1[ca,ch,production], -g1_mix*gi_mix/gi_ggHBSM)]
                                                 for ca, ch, production in itertools.product(categories, channels, config.productionsforcombine)),
                                               []
                                             ))

    VVHSM     = self.DbkgSum("VVHSM", *sum(
                                           ([(allVBFg14gi0[ca,ch,production], 1), (allZHg14gi0[ca,ch,production], 1), (allWHg14gi0[ca,ch,production], 1)]
                                              for ca, ch, production in itertools.product(categories, channels, config.productionsforcombine)),
                                            []
                                          ))
    VVHBSM    = self.DbkgSum("VVHBSM", *sum(
                                            ([(allVBFg10gi4[ca,ch,production], 1), (allZHg10gi4[ca,ch,production], 1), (allWHg10gi4[ca,ch,production], 1)]
                                               for ca, ch, production in itertools.product(categories, channels, config.productionsforcombine)),
                                             []
                                           ))
    VVHmix_p  = self.DbkgSum("VVHmix_p", *sum((
                                                [(template, g1_mix**(4-j) * (+gi_mix)**j) for j, template in enumerate(allVBFpieces[ca,ch,production])]
                                              + [(template, g1_mix**(4-j) * (+gi_mix)**j) for j, template in enumerate(allZHpieces[ca,ch,production])]
                                              + [(template, g1_mix**(4-j) * (+gi_mix)**j) for j, template in enumerate(allWHpieces[ca,ch,production])]
                                                 for ca, ch, production in itertools.product(categories, channels, config.productionsforcombine)),
                                               []
                                             ))
    VVHmix_m  = self.DbkgSum("VVHmix_m", *sum((
                                                [(template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(allVBFpieces[ca,ch,production])]
                                              + [(template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(allZHpieces[ca,ch,production])]
                                              + [(template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(allWHpieces[ca,ch,production])]
                                                 for ca, ch, production in itertools.product(categories, channels, config.productionsforcombine)),
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

    for ca, ch, production in itertools.product(categories, channels, config.productionsforcombine):
      ggHg12gi0[ca,ch,production] = self.TemplateFromFile(   "ggH", ca, self.enrichstatus, self.normalization, production, ch, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
      ggHg10gi2[ca,ch,production] = self.TemplateFromFile(   "ggH", ca, self.enrichstatus, self.normalization, production, ch, self.shapesystematic, self.analysis, BSMhypothesis)
      ggHg11gi1[ca,ch,production] = self.IntTemplateFromFile("ggH", ca, self.enrichstatus, self.normalization, production, ch, self.shapesystematic, self.analysis, "g11gi1")

      ttHg12gi0[ca,ch,production] = self.TemplateFromFile(   "ttH", ca, self.enrichstatus, self.normalization, production, ch, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
      ttHg10gi2[ca,ch,production] = self.TemplateFromFile(   "ttH", ca, self.enrichstatus, self.normalization, production, ch, self.shapesystematic, self.analysis, BSMhypothesis)
      ttHg11gi1[ca,ch,production] = self.IntTemplateFromFile("ttH", ca, self.enrichstatus, self.normalization, production, ch, self.shapesystematic, self.analysis, "g11gi1")

      VBFg14gi0[ca,ch,production] = self.TemplateFromFile(   "VBF", ca, self.enrichstatus, self.normalization, production, ch, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
      VBFg13gi1[ca,ch,production] = self.IntTemplateFromFile("VBF", ca, self.enrichstatus, self.normalization, production, ch, self.shapesystematic, self.analysis, "g13gi1")
      VBFg12gi2[ca,ch,production] = self.IntTemplateFromFile("VBF", ca, self.enrichstatus, self.normalization, production, ch, self.shapesystematic, self.analysis, "g12gi2")
      VBFg11gi3[ca,ch,production] = self.IntTemplateFromFile("VBF", ca, self.enrichstatus, self.normalization, production, ch, self.shapesystematic, self.analysis, "g11gi3")
      VBFg10gi4[ca,ch,production] = self.TemplateFromFile(   "VBF", ca, self.enrichstatus, self.normalization, production, ch, self.shapesystematic, self.analysis, BSMhypothesis)

      VBFpieces[ca,ch,production] = [VBFg14gi0[ca,ch,production], VBFg13gi1[ca,ch,production], VBFg12gi2[ca,ch,production], VBFg11gi3[ca,ch,production], VBFg10gi4[ca,ch,production]]

      ZHg14gi0[ca,ch,production] = self.TemplateFromFile(   "ZH", ca, self.enrichstatus, self.normalization, production, ch, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
      ZHg13gi1[ca,ch,production] = self.IntTemplateFromFile("ZH", ca, self.enrichstatus, self.normalization, production, ch, self.shapesystematic, self.analysis, "g13gi1")
      ZHg12gi2[ca,ch,production] = self.IntTemplateFromFile("ZH", ca, self.enrichstatus, self.normalization, production, ch, self.shapesystematic, self.analysis, "g12gi2")
      ZHg11gi3[ca,ch,production] = self.IntTemplateFromFile("ZH", ca, self.enrichstatus, self.normalization, production, ch, self.shapesystematic, self.analysis, "g11gi3")
      ZHg10gi4[ca,ch,production] = self.TemplateFromFile(   "ZH", ca, self.enrichstatus, self.normalization, production, ch, self.shapesystematic, self.analysis, BSMhypothesis)

      ZHpieces[ca,ch,production] = [ZHg14gi0[ca,ch,production], ZHg13gi1[ca,ch,production], ZHg12gi2[ca,ch,production], ZHg11gi3[ca,ch,production], ZHg10gi4[ca,ch,production]]

      WHg14gi0[ca,ch,production] = self.TemplateFromFile(   "WH", ca, self.enrichstatus, self.normalization, production, ch, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
      WHg13gi1[ca,ch,production] = self.IntTemplateFromFile("WH", ca, self.enrichstatus, self.normalization, production, ch, self.shapesystematic, self.analysis, "g13gi1")
      WHg12gi2[ca,ch,production] = self.IntTemplateFromFile("WH", ca, self.enrichstatus, self.normalization, production, ch, self.shapesystematic, self.analysis, "g12gi2")
      WHg11gi3[ca,ch,production] = self.IntTemplateFromFile("WH", ca, self.enrichstatus, self.normalization, production, ch, self.shapesystematic, self.analysis, "g11gi3")
      WHg10gi4[ca,ch,production] = self.TemplateFromFile(   "WH", ca, self.enrichstatus, self.normalization, production, ch, self.shapesystematic, self.analysis, BSMhypothesis)

      WHpieces[ca,ch,production] = [WHg14gi0[ca,ch,production], WHg13gi1[ca,ch,production], WHg12gi2[ca,ch,production], WHg11gi3[ca,ch,production], WHg10gi4[ca,ch,production]]

    del ca, ch, production

#    ggHSM = self.TemplateSum("ggH SM", (ggHg12gi0[category,channel,production], 1))
#    ggHBSM = self.ComponentTemplateSumInGroup("ggH {}=1".format(fainame), ffHBSM, ffHSM, (ggHg10gi2[category,channel,production], 1))
#    ggHmix_p  = self.ComponentTemplateSumInGroup("ggH {}=#plus0.5" .format(fainame), ffHmix_p, ffHSM,
#                                                 (ggHg12gi0[category,channel], g1_mix**2), (ggHg10gi2[category,channel], (gi_mix/gi_ggHBSM)**2), (ggHg11gi1[category,channel],  g1_mix*gi_mix/gi_ggHBSM))
#    ggHmix_m  = self.ComponentTemplateSumInGroup("ggH {}=#minus0.5".format(fainame), ffHmix_m, ffHSM,
#                                                 (ggHg12gi0[category,channel], g1_mix**2), (ggHg10gi2[category,channel], (gi_mix/gi_ggHBSM)**2), (ggHg11gi1[category,channel], -g1_mix*gi_mix/gi_ggHBSM))
#    ttHSM = self.TemplateSum("ttH SM", (ttHg12gi0[category,channel], 1))
#    ttHBSM = self.ComponentTemplateSumInGroup("ttH {}=1".format(fainame), ffHBSM, ffHSM, (ttHg10gi2[category,channel], 1))
#    ttHmix_p  = self.ComponentTemplateSumInGroup("ttH {}=#plus0.5" .format(fainame), ffHmix_p, ffHSM,
#                                                 (ttHg12gi0[category,channel], g1_mix**2), (ttHg10gi2[category,channel], (gi_mix/gi_ggHBSM)**2), (ttHg11gi1[category,channel],  g1_mix*gi_mix/gi_ggHBSM))
#    ttHmix_m  = self.ComponentTemplateSumInGroup("ttH {}=#minus0.5".format(fainame), ffHmix_m, ffHSM,
#                                                 (ttHg12gi0[category,channel], g1_mix**2), (ttHg10gi2[category,channel], (gi_mix/gi_ggHBSM)**2), (ttHg11gi1[category,channel], -g1_mix*gi_mix/gi_ggHBSM))

#    VBFSM = self.TemplateSum("VBF SM", (VBFg14gi0[category,channel], 1))
#    VBFBSM = self.ComponentTemplateSumInGroup("VBF {}=1".format(fainame), VVHBSM, VVHSM, (VBFg10gi4[category,channel], 1))
#    VBFmix_p  = self.ComponentTemplateSumInGroup("VBF {}=#plus0.5" .format(fainame), VVHmix_p, VVHSM,
#                                            *((template, g1_mix**(4-j) * gi_mix**j) for j, template in enumerate(VBFpieces[category,channel]))
#                                           )
#    VBFmix_m  = self.ComponentTemplateSumInGroup("VBF {}=#minus0.5".format(fainame), VVHmix_m, VVHSM,
#                                            *((template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(VBFpieces[category,channel]))
#                                           )

#    VHSM = self.TemplateSum("VH SM", (ZHg14gi0[category,channel], 1), (WHg14gi0[category,channel], 1))
#    VHBSM = self.ComponentTemplateSumInGroup("VH {}=1".format(fainame), VVHBSM, VVHSM, (ZHg10gi4[category,channel], 1), (WHg10gi4[category,channel], 1))
#    VHmix_p  = self.ComponentTemplateSumInGroup("VH {}=#plus0.5" .format(fainame), VVHmix_p, VVHSM,
#                                            *(
#                                                [(template, g1_mix**(4-j) * gi_mix**j) for j, template in enumerate(ZHpieces[category,channel])]
#                                              + [(template, g1_mix**(4-j) * gi_mix**j) for j, template in enumerate(WHpieces[category,channel])]
#                                             )
#                                           )
#    VHmix_m  = self.ComponentTemplateSumInGroup("VH {}=#minus0.5" .format(fainame), VVHmix_m, VVHSM,
#                                            *(
#                                                [(template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(ZHpieces[category,channel])]
#                                              + [(template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(WHpieces[category,channel])]
#                                             )
#                                           )

#    SM    = self.TemplateSum("SM",                                 (ggHSM,    ggHfactor), (VBFSM,    VBFfactor), (VHSM,    VHfactor), (ttHSM,    ttHfactor), linecolor=1)
#    BSM   = self.TemplateSum("{}=1".format(self.analysis.title()), (ggHBSM,   ggHfactor), (VBFBSM,   VBFfactor), (VHBSM,   VHfactor), (ttHBSM,   ttHfactor), linecolor=ROOT.kCyan)
#    mix_p = self.TemplateSum("{}=#plus0.5".format(fainame),        (ggHmix_p, ggHfactor), (VBFmix_p, VBFfactor), (VHmix_p, VHfactor), (ttHmix_p, ttHfactor), linecolor=ROOT.kGreen+3)
#    mix_m = self.TemplateSum("{}=#minus0.5".format(fainame),       (ggHmix_m, ggHfactor), (VBFmix_m, VBFfactor), (VHmix_m, VHfactor), (ttHmix_m, ttHfactor), linecolor=4)

#    qqZZ      = {production: self.TemplateFromFile(category, self.enrichstatus, self.normalization, self.analysis, channel, "qqZZ",    self.shapesystematic, production, linecolor=6) for production in config.productionsforcombine}
#    ggZZ      = {production: self.TemplateFromFile(category, self.enrichstatus, self.normalization, self.analysis, channel, "ggZZ",    self.shapesystematic, production, linecolor=ROOT.kOrange+6) for production in config.productionsforcombine}
#    VBFbkg    = {production: self.TemplateFromFile(category, self.enrichstatus, self.normalization, self.analysis, channel, "VBF bkg", self.shapesystematic, production, linecolor=ROOT.kViolet+7) for production in config.productionsforcombine}
#    ZX        = {production: self.TemplateFromFile(category, self.enrichstatus, self.normalization, self.analysis, channel, "ZX",      self.shapesystematic, production, linecolor=2) for production in config.productionsforcombine}

    templates = []
    if not animation and not nicestyle:
        templates += [
                      SM, BSM, mix_p, mix_m,
                     ]
    if animation and not nicestyle:
        assert ggHfactor == VBFfactor == VHfactor == 1
        if customfai <  0: plusminus = "#minus"
        if customfai == 0: plusminus = ""
        if customfai  > 0: plusminus = "#plus"
        ggHcustom = self.ComponentTemplateSum("ggH ({}={}{:.2f})".format(self.analysis.title(superscript="dec"), plusminus, abs(customsample.fai("ggH", self.analysis))), 1, (ggHg12gi0[category,channel], g1_custom**2), (ggHg10gi2[category,channel], (gi_custom/gi_ggHBSM)**2), (ggHg11gi1[category,channel],  g1_custom*gi_custom/gi_ggHBSM), linecolor=1)
        VBFcustom = self.ComponentTemplateSum("VBF ({}={}{:.2f})".format(self.analysis.title(superscript="VBF"), plusminus, abs(customsample.fai("VBF", self.analysis))), 1,
                                              *((template, g1_custom**(4-j) * gi_custom**j) for j, template in enumerate(VBFpieces[category,channel])),
                                              linecolor=2
                                             )
        VHcustom  = self.ComponentTemplateSum("VH ({}={}{:.2f})".format(self.analysis.title(superscript="VH"), plusminus, abs(customsample.fai("VH", self.analysis))), 1,
                                              *(
                                                  [(template, g1_custom**(4-j) * gi_custom**j) for j, template in enumerate(ZHpieces[category,channel])]
                                                + [(template, g1_custom**(4-j) * gi_custom**j) for j, template in enumerate(WHpieces[category,channel])]
                                               ),
                                               linecolor=4
                                             )
        for t in ggHcustom, VBFcustom, VHcustom:
            assert abs(t.h.Integral() - 1) < 1e-6, t.h.Integral()

        templates += [
                      ggHcustom, VBFcustom, VHcustom
                     ]

    if not justoneproductionmode and not animation and not nicestyle:
        templates += [
                      qqZZ, ggZZ, VBFbkg, ZX,
                     ]

    if nicestyle:
        ZXname = "Z+X"
        ZXkwargs = dict(linecolor=1, fillcolor=ROOT.TColor.GetColor("#669966"), linewidth=2, fillstyle=1001, legendoption="f")
        ZZname = "ZZ/Z#gamma*"
        ZZkwargs = dict(linecolor=1, fillcolor=ROOT.kAzure-9, linewidth=2, fillstyle=1001, legendoption="f")
        if Dbkg_allcategories:
            ZX = self.DbkgSum(ZXname,
                              *((self.TemplateFromFile(ca, self.enrichstatus, self.normalization, self.analysis, ch, "ZX",      self.shapesystematic, production, linecolor=2), 1) for ca, ch, production in itertools.product(categories, channels, config.productionsforcombine)),
                              **ZXkwargs)
            ZZ = self.DbkgSum(ZZname,
                              *((self.TemplateFromFile(ca, self.enrichstatus, self.normalization, self.analysis, ch, p,      self.shapesystematic, production, linecolor=2), 1) for ca, ch, production in itertools.product(categories, channels, config.productionsforcombine) for p in ("qqZZ", "ggZZ", "VBFbkg", "ZX")),
                              **ZZkwargs)
        else:
            ZX = self.TemplateSum(ZXname,
                                  *((self.TemplateFromFile(category, self.enrichstatus, self.normalization, self.analysis, ch, "ZX",      self.shapesystematic, production, linecolor=2), 1) for ch, production in itertools.product(channels, config.productionsforcombine)),
                                  **ZXkwargs)
            ZZ = self.TemplateSum(ZZname,
                                  *((self.TemplateFromFile(category, self.enrichstatus, self.normalization, self.analysis, ch, p,      self.shapesystematic, production, linecolor=2), 1) for ch in channels for production in config.productionsforcombine for p in ("qqZZ", "ggZZ", "VBFbkg", "ZX")),
                                  **ZZkwargs)
        if category == "Untagged":
            superscript = None
        elif category == "VBFtagged":
            superscript = "VBF"
        elif category == "VHHadrtagged":
            superscript = "VH"

        templates += [ZZ, ZX]

        legendargs = [(0.20,0.57,0.58,0.90), (0.20,0.57,0.58,0.90), (0.23,0.57,0.61,0.90)]
        if category == "VBFtagged" and self.analysis == "fa2":
            legendargs[2] = (0.26,0.57,0.64,0.90)
        if Dbkg_allcategories:
            usecategories = categories
            templateclass = self.DbkgSum
            componenttemplateclass = self.ComponentDbkgSum
            componenttemplateclassingroup = self.ComponentDbkgSumInGroup
        else:
            usecategories = [category]
            templateclass = self.TemplateSum
            componenttemplateclass = self.ComponentTemplateSum
            componenttemplateclassingroup = self.ComponentTemplateSumInGroup

        if not animation:
            SMffH = templateclass("",
                                  *sum(
                                       ([(ggHg12gi0[ca,ch,production], 1), (ttHg12gi0[ca,ch,production], 1)]
                                           for ca, ch, production in itertools.product(usecategories, channels, config.productionsforcombine)),
                                        []
                                      )
                                 )
            SMVVH = templateclass("",
                                  *sum(
                                       ([(VBFg14gi0[ca,ch,production], 1), (ZHg14gi0[ca,ch,production], 1), (WHg14gi0[ca,ch,production], 1)]
                                           for ca, ch, production in itertools.product(usecategories, channels, config.productionsforcombine)),
                                        []
                                      )
                                 )
            BSMffH = componenttemplateclass("", SMffH,
                                          *sum(
                                               ([(ggHg10gi2[ca,ch,production], 1), (ttHg10gi2[ca,ch,production], 1)]
                                                   for ca, ch, production in itertools.product(usecategories, channels, config.productionsforcombine)),
                                                []
                                               )
                                         )
            BSMVVH = componenttemplateclass("", SMVVH,
                                          *sum(
                                               ([(VBFg10gi4[ca,ch,production], 1), (ZHg10gi4[ca,ch,production], 1), (WHg10gi4[ca,ch,production], 1)]
                                                   for ca, ch, production in itertools.product(usecategories, channels, config.productionsforcombine)),
                                                []
                                               )
                                         )
            mix_pffH = componenttemplateclass("", SMffH,
                                            *sum(
                                                 ([(ggHg12gi0[ca,ch,production], g1_mix**2), (ggHg10gi2[ca,ch,production], (gi_mix/gi_ggHBSM)**2), (ggHg11gi1[ca,ch,production],  g1_mix*gi_mix/gi_ggHBSM),
                                                   (ttHg12gi0[ca,ch,production], g1_mix**2), (ttHg10gi2[ca,ch,production], (gi_mix/gi_ggHBSM)**2), (ttHg11gi1[ca,ch,production],  g1_mix*gi_mix/gi_ggHBSM)]
                                                    for ca, ch, production in itertools.product(usecategories, channels, config.productionsforcombine)),
                                                  []
                                                 )
                                           )
            mix_pVVH = componenttemplateclass("", SMVVH,
                                            *sum(
                                                 ([(template, g1_mix**(4-j) * (+gi_mix)**j) for j, template in enumerate(VBFpieces[ca,ch,production])]
                                                + [(template, g1_mix**(4-j) * (+gi_mix)**j) for j, template in enumerate(ZHpieces[ca,ch,production])]
                                                + [(template, g1_mix**(4-j) * (+gi_mix)**j) for j, template in enumerate(WHpieces[ca,ch,production])]
                                                    for ca, ch, production in itertools.product(usecategories, channels, config.productionsforcombine)),
                                                  []
                                                 )
                                           )
            mix_mffH = componenttemplateclass("", SMffH,
                                            *sum(
                                                 ([(ggHg12gi0[ca,ch,production], g1_mix**2), (ggHg10gi2[ca,ch,production], (gi_mix/gi_ggHBSM)**2), (ggHg11gi1[ca,ch,production], -g1_mix*gi_mix/gi_ggHBSM),
                                                   (ttHg12gi0[ca,ch,production], g1_mix**2), (ttHg10gi2[ca,ch,production], (gi_mix/gi_ggHBSM)**2), (ttHg11gi1[ca,ch,production], -g1_mix*gi_mix/gi_ggHBSM)]
                                                    for ca, ch, production in itertools.product(usecategories, channels, config.productionsforcombine)),
                                                  []
                                                 )
                                           )
            mix_mVVH = componenttemplateclass("", SMVVH,
                                            *sum(
                                                 ([(template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(VBFpieces[ca,ch,production])]
                                                + [(template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(ZHpieces[ca,ch,production])]
                                                + [(template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(WHpieces[ca,ch,production])]
                                                    for ca, ch, production in itertools.product(usecategories, channels, config.productionsforcombine)),
                                                  []
                                                 )
                                           )

            SMbottom = SMVVH
            BSMbottom = BSMVVH
            mix_pbottom = mix_pVVH
            mix_mbottom = mix_mVVH
            bottomtitle = "VBF+VH"
            bottomcolor = 4
            bottommu = muV
            #toptitle = "ggH+t#bar{t}H"
            toptitle = "Total"
            topcolor = ROOT.kOrange+10

            SMbottom = self.TemplateSum("{} SM".format(bottomtitle),
                                       (SMbottom, bottommu), (ZZ, 1), #which already has ZX
                                       linecolor=bottomcolor, linewidth=2, legendoption="f")
            BSMbottom = self.TemplateSum("{} {} = 1".format(bottomtitle, self.analysis.title()),
                                         (BSMbottom, bottommu), (ZZ, 1),
                                         linecolor=bottomcolor, linewidth=2, linestyle=2, legendoption="f")
            mix_pbottom = self.TemplateSum("{} {} = #plus 0.5".format(bottomtitle, self.analysis.title(superscript=superscript)),
                                           (mix_pbottom, bottommu), (ZZ, 1),
                                           linecolor=bottomcolor, linewidth=2, linestyle=2, legendoption="f")
            mix_mbottom = self.TemplateSum("{} {} = #minus 0.5".format(bottomtitle, self.analysis.title(superscript=superscript)),
                                           (mix_mbottom, bottommu), (ZZ, 1),
                                           linecolor=bottomcolor, linewidth=2, linestyle=2, legendoption="f")

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
            templates[0:0] = [SM, SMbottom, BSM, BSMbottom, mix_p, mix_pbottom, mix_m, mix_mbottom] #will remove some later, depending on the discriminant

            rebin = 1

        else: #animation
            ffHintegral  = self.DbkgSum("ffHintegral", *sum(
                                                            ([(allggHg12gi0[ca,ch], g1_custom**2), (allggHg10gi2[ca,ch], (gi_custom/gi_ggHBSM)**2), (allggHg11gi1[ca,ch],  g1_custom*gi_custom/gi_ggHBSM),
                                                            (allttHg12gi0[ca,ch], g1_custom**2), (allttHg10gi2[ca,ch], (gi_custom/gi_ggHBSM)**2), (allttHg11gi1[ca,ch],  g1_custom*gi_custom/gi_ggHBSM)]
                                                             for ca, ch, production in itertools.product(categories, channels, config.productionsforcombine)),
                                                            []
                                                           ))
            VVHintegral  = self.DbkgSum("VVHintegral", *sum((
                                                              [(template, g1_custom**(4-j) * (+gi_custom)**j) for j, template in enumerate(allVBFpieces[ca,ch])]
                                                            + [(template, g1_custom**(4-j) * (+gi_custom)**j) for j, template in enumerate(allZHpieces[ca,ch])]
                                                            + [(template, g1_custom**(4-j) * (+gi_custom)**j) for j, template in enumerate(allWHpieces[ca,ch])]
                                                               for ca, ch, production in itertools.product(categories, channels, config.productionsforcombine)),
                                                              []
                                                           ))
            ffH = componenttemplateclassingroup("", ffHintegral, ffHSM,
                                       *sum(
                                            ([(ggHg12gi0[ca,ch], g1_custom**2), (ggHg10gi2[ca,ch], (gi_custom/gi_ggHBSM)**2), (ggHg11gi1[ca,ch],  g1_custom*gi_custom/gi_ggHBSM),
                                              (ttHg12gi0[ca,ch], g1_custom**2), (ttHg10gi2[ca,ch], (gi_custom/gi_ggHBSM)**2), (ttHg11gi1[ca,ch],  g1_custom*gi_custom/gi_ggHBSM)]
                                               for ca, ch in itertools.product(usecategories, channels)),
                                             []
                                            )
                                      )
            VVH = componenttemplateclassingroup("", VVHintegral, VVHSM,
                                       *sum(
                                            ([(template, g1_custom**(4-j) * (+gi_custom)**j) for j, template in enumerate(VBFpieces[ca,ch])]
                                           + [(template, g1_custom**(4-j) * (+gi_custom)**j) for j, template in enumerate(ZHpieces[ca,ch])]
                                           + [(template, g1_custom**(4-j) * (+gi_custom)**j) for j, template in enumerate(WHpieces[ca,ch])]
                                               for ca, ch in itertools.product(usecategories, channels)),
                                             []
                                            )
                                      )

            bottom = VVH
            bottomtitle = "VBF+VH"
            bottomcolor = 4
            bottommu = muV
            toptitle = "Total" #"ggH+t#bar{t}H"
            topcolor = ROOT.kOrange+10

            bottom = self.TemplateSum(bottomtitle,
                                     (bottom, bottommu), (ZZ, 1), #which already has ZX
                                     linecolor=bottomcolor, linewidth=2, legendoption="f")
            top = self.TemplateSum(toptitle,
                                   (VVH, muV), (ffH, muf), (ZZ, 1), #which already has ZX
                                   linecolor=topcolor, linewidth=2)

            templates[0:0] = [top, bottom] #will remove some later, depending on the discriminant

            rebin = 1

        if self.enrichstatus == "impoverish" and config.showblinddistributions or config.unblinddistributions:
            if Dbkg_allcategories:
                data = self.DbkgSum("",
                                    *((self.TemplateFromFile(
                                                             ca, self.enrichstatus, self.normalization,
                                                             self.analysis, ch, "data", self.production,
                                                             linecolor=1,
                                                            ), 1) for ca, ch, production in itertools.product(categories, channels, config.productionsforcombine))
                                   )
                data = [None, None, style.asymmerrorsfromhistogram(data.Projection(2, rebin=rebin), showemptyerrors=False)]
            else:
                data = self.TemplateSum("",
                                        *((self.TemplateFromFile(
                                                                 category, self.enrichstatus, self.normalization,
                                                                 self.analysis, ch, "data", production,
                                                                 linecolor=1,
                                                                ), 1) for ch in channels for production in config.productionsforcombine)
                                       )
                data = [style.asymmerrorsfromhistogram(data.Projection(i, rebin=rebin), showemptyerrors=False) for i in range(3)]
            for g in data:
                if g is None: continue
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
        if Dbkg_allcategories and i != 2: continue
        usetemplates = templates
        if nicestyle and not animation:
            usetemplates = templates[:]
            tf = TemplatesFile("2e2mu", "ggh", category, config.productionsforcombine[0], self.analysis)
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
            if category == "Untagged" and discriminant.name == "D_bkg" and self.enrichstatus == "fullrange":
                hstack.SetMaximum(ymax)
#            elif category == "VBFtagged" and discriminant.name == "D_bkg" and self.enrichstatus == "fullrange":
#                if self.analysis in ("fa2", "fL1", "fL1Zg"):
#                    hstack.SetMaximum(ymax * 1.25)
#                else:
#                    hstack.SetMaximum(ymax)
            elif category == "VHHadrtagged" and discriminant.name == "D_bkg" and self.enrichstatus == "fullrange":
                if animation:
                    hstack.SetMaximum(ymax * 1.7)
                else:
                    hstack.SetMaximum(ymax)
            elif discriminant.name in ("D_0minus_decay", "D_CP_decay_new", "D_0minus_HadVHdecay", "D_0hplus_decay") and self.enrichstatus == "enrich":
                hstack.SetMaximum(ymax * 1.4)
            elif discriminant.name in ("D_L1_decay", "D_L1Zg_decay") and self.enrichstatus == "enrich" and not animation:
                hstack.SetMaximum(ymax * 1.4)
            elif discriminant.name in ("D_CP_VBF_new", "D_0hplus_VBFdecay", "D_int_VBF_new", "D_int_HadVH_new") and self.enrichstatus == "enrich":
                hstack.SetMaximum(ymax * 1.7)
            elif discriminant.name == "D_CP_HadVH_new" and self.enrichstatus == "enrich":
                hstack.SetMaximum(ymax * 1.7)
            elif discriminant.name == "D_0hplus_HadVHdecay" and self.enrichstatus == "enrich":
                if animation:
                    hstack.SetMaximum(ymax * 1.7)
                else:
                    hstack.SetMaximum(ymax * 1.4)
            elif discriminant.name in ("D_L1_VBFdecay", "D_L1_HadVHdecay", "D_L1Zg_VBFdecay", "D_L1Zg_HadVHdecay") and self.enrichstatus == "enrich" and animation:
                hstack.SetMaximum(ymax * 1.7)
            else:
                hstack.SetMaximum(ymax * 1.25)

        c1.cd()
        hstack.Draw("nostack")
        hstack.GetXaxis().SetTitle(discriminant.title)
        hstack.GetYaxis().SetTitle(
#                                   "Events / {:.2f}".format(
#                                                            (hstack.GetXaxis().GetXmax() - hstack.GetXaxis().GetXmin()) / hstack.GetXaxis().GetNbins()
#                                                           )
                                   "Events / bin"
                                  )
        style.applyaxesstyle(hstack)
        if nicestyle:
            subfigletter = None
            #aux - add subfig letters here?
            if forWIN:
                CMStext = "Preliminary"
            elif discriminant.name == "D_bkg_VBFdecay" and not Dbkg_allcategories and self.enrichstatus == "fullrange" and self.analysis == "fa3" and not animation:
                CMStext = ""
            elif discriminant.name == "D_bkg_HadVHdecay" and not Dbkg_allcategories and self.enrichstatus == "fullrange" and self.analysis == "fa3" and not animation:
                CMStext = ""
            elif discriminant.name == "D_bkg" and not Dbkg_allcategories and self.enrichstatus == "fullrange" and self.analysis == "fa3" and not animation:
                CMStext = ""
            elif discriminant.name == "D_0minus_decay" and self.enrichstatus == "enrich" and not animation:
                CMStext = ""
            elif discriminant.name == "D_0hplus_decay" and self.enrichstatus == "enrich" and not animation and self.analysis == "fa2":
                CMStext = ""
            elif discriminant.name == "D_L1_decay" and self.enrichstatus == "enrich" and not animation:
                CMStext = ""
            elif discriminant.name == "D_L1Zg_decay" and self.enrichstatus == "enrich" and not animation:
                CMStext = ""
            elif discriminant.name == "D_0minus_VBFdecay" and self.enrichstatus == "enrich" and not animation:
                CMStext = ""
            elif discriminant.name == "D_0minus_HadVHdecay" and self.enrichstatus == "enrich" and not animation:
                CMStext = ""
            elif discriminant.name == "D_CP_decay_new" and self.enrichstatus == "enrich" and not animation:
                CMStext = ""
            elif discriminant.name in ("D_bkg", "D_bkg_VBFdecay", "D_bkg_HadVHdecay") and not Dbkg_allcategories and self.enrichstatus == "fullrange":
                CMStext = "Supplementary"
            elif discriminant.name not in ("D_bkg", "D_bkg_VBFdecay", "D_bkg_HadVHdecay") and self.enrichstatus == "enrich":
                CMStext = "Supplementary"
            else:
                CMStext = "Internal"

            lumi = sum(float(Luminosity("fordata", production)) for production in config.productionsforcombine)

            style.CMS(CMStext, lumi)

            hstack.GetXaxis().CenterTitle()
            hstack.GetYaxis().CenterTitle()

            cuttextkwargs = {}
            if animation:
                cuttextkwargs.update(x1=.48+.03*(discriminant.name=="D_bkg"), x2=.58+.03*(discriminant.name=="D_bkg"))
            if subfigletter is not None:
                style.subfig(subfigletter)
                cuttextkwargs.update(y1=.78, y2=.86)
            style.cuttext(category.nameonplot, **cuttextkwargs)

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

  def saveasdir(self, *categoryandchannel, **kwargs):
      forWIN = kwargs.pop("forWIN", False)
      assert not kwargs
      forWIN = "forPAS" if forWIN else ""
      categoryandchannel = self.CategoryAndChannel(*categoryandchannel)
      assert self.normalization == "rescalemixtures"
      return os.path.join(config.plotsbasedir, "templateprojections", forWIN, "projections", self.enrichstatus.dirname(), "{}_{}/{}/{}".format(self.analysis, self.production, categoryandchannel.category, categoryandchannel.channel))
  def saveasdir_niceplots(self, category, forWIN=False):
      assert self.normalization == "rescalemixtures"# and len(config.productionsforcombine) == 1
      forWIN = "forPAS" if forWIN else ""
      result = os.path.join(config.plotsbasedir, "templateprojections", forWIN, "niceplots", self.enrichstatus.dirname(), "{}/{}".format(self.analysis, Category(category)))
      return result
  def saveasdir_Dbkgsum(self, forWIN=False):
      forWIN = "forPAS" if forWIN else ""
      return os.path.join(config.plotsbasedir, "templateprojections", forWIN, "niceplots", self.enrichstatus.dirname(), str(self.analysis))

  def discriminants(self, category):
      return TemplatesFile("2e2mu", self.shapesystematic, "ggh", self.analysis, config.productionsforcombine[0], category).discriminants

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
        pt.AddText("{}={:.2f}".format(self.analysis.title(), self.fai_decay))
        pt.AddText("{}={:.2f}".format(self.analysis.title(superscript="VBF"), self.fai("VBF", self.analysis)))
        pt.AddText("{}={:.2f}".format(self.analysis.title(superscript="VH"), self.fai("VH", self.analysis)))
#        pt.AddText("{}={:.2f}".format("#mu_{V}", self.muV))
#        pt.AddText("{}={:.2f}".format("#mu_{f}", self.muf))
        pt.AddText("{}={:.2f}".format("-2#Deltaln L", self.deltaNLL))
        return pt

  def animation(self, category, channel, floor=False, nicestyle=False, Dbkg_allcategories=False, forWIN=False):
    with mkdtemp() as tmpdir:
      category = Category(category)
      channel = Channel(channel)
  
      nsteps = 200
  
      finaldir = os.path.join(self.saveasdir(category, channel, forWIN=forWIN), "animation")
  
      if nicestyle:
        if channel != "2e2mu": return
        if Dbkg_allcategories:
          if category != "Untagged": return
          finaldir = self.saveasdir_Dbkgsum(forWIN=forWIN)
        else:
          finaldir = os.path.join(self.saveasdir_niceplots(category, forWIN=forWIN), "animation")
        animation = self.animationstepsforniceplots(self.analysis)
  
      else:
        assert not Dbkg_allcategories
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
        "nicestyle": nicestyle,
        "Dbkg_allcategories": Dbkg_allcategories,
        "forWIN": forWIN,
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
  
      for j, discriminant in enumerate(self.discriminants(category)):
        if Dbkg_allcategories and j != 2: continue
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
        finalplot = os.path.join(finaldir, "{}.gif".format(discriminant.name))
        convertcommand.append(finalplot)
        #http://stackoverflow.com/a/38792806/5228524
        #subprocess.check_call(convertcommand)
        os.system(" ".join(pipes.quote(_) for _ in convertcommand))

  @classmethod
  def scantreeforanimations(cls, analysis):
    analysis = Analysis(analysis)
    scantree = ROOT.TChain("limit")
    for filename in glob.glob(os.path.join(config.repositorydir, "CMSSW_8_1_0", "src", "HiggsAnalysis", "HZZ4l_Combination",
                                       "CreateDatacards", "cards_{}_allsysts".format(analysis), "higgsCombine_obs_*.root")):
        if re.match("higgsCombine_obs_lumi[0-9.]*_7813(_[0-9]*,[-.0-9]*,[-.0-9]*)*.MultiDimFit.mH125.root", os.path.basename(filename)):
            scantree.Add(filename)
    return scantree

  @classmethod
  def animationstepsforniceplots(cls, analysis):
    animation = []
    scantree = cls.scantreeforanimations(analysis)
    minNLL, faiforminNLL = min((scantree.deltaNLL+scantree.nll+scantree.nll0, scantree.CMS_zz4l_fai1) for entry in scantree)
    VBFmix = samplewithfai("ggH", analysis, 0.5, productionmodeforfai="VBF").fai("ggH", analysis)
    VHmix = samplewithfai("ggH", analysis, 0.5, productionmodeforfai="VH").fai("ggH", analysis)
    pauses = [
              min((scantree.CMS_zz4l_fai1 for entry in scantree), key = lambda x: abs(x-_))
                 for _ in (0, -1, 1, .5, -.5, VBFmix, -VBFmix, VHmix, -VHmix) if _ in (0, 1, -1)
             ]
    for entry in scantree:
        animation.append(
                         cls.AnimationStep(
                                           "ggH",
                                           analysis,
                                           scantree.CMS_zz4l_fai1,
                                           100 if scantree.CMS_zz4l_fai1 == faiforminNLL or scantree.CMS_zz4l_fai1 in pauses
                                              else 10,
                                           muV = scantree.muV_scaled,
                                           muf = scantree.muf_scaled,
                                           deltaNLL = 2*(scantree.deltaNLL+scantree.nll+scantree.nll0 - minNLL)
                                          )
                        )
    animation = sorted(set(animation))
    while True:
        for step, nextstep in pairwise(animation[:]):
            if step.fai_decay == nextstep.fai_decay:
                animation.remove(step)
                break
        else:
            break
    return animation




def projections(*args):
    Projections(*args).projections()

def main():
  if sys.argv[1] == "submitjobs":
    if sys.argv[2:]:
      raise ValueError("For submitjobs, don't add any more arguments")
    for process in 1, 2, 4, 5:
      for analysis in analyses:
        for enrichstatus in enrichstatuses:
          if enrichstatus == "blind": continue

          words = [
            "time",
            os.path.join(config.repositorydir, "step10_plottingutilities", "projections.py"),
            str(process),
            str(analysis),
          ]

          if process in (4, 5):
            words.append(str(enrichstatus))
          elif enrichstatus != "enrich":
            continue

          jobtext = " ".join(pipes.quote(_) for _ in words)

          submitjobkwargs = {
            "jobname": "{} {}".format(process, analysis),
            "email": True,
          }
          if process in (1, 2):
            submitjobkwargs["jobtime"] = "0:20:0"
          elif process in (4, 5):
            if config.host == "MARCC":
              submitjobkwargs["queue"] = "lrgmem"
              submitjobkwargs["memory"] = "10000M"
            submitjobkwargs["jobtime"] = "1:0:0"

          submitjob(jobtext, **submitjobkwargs)

    return

  useanalyses = []
  useenrichstatuses = []
  for _ in sys.argv[2:]:
    try:
      useanalyses.append(Analysis(_))
    except ValueError:
      useenrichstatuses.append(EnrichStatus(_))
  if not useanalyses:
    useanalyses = analyses
  if not useenrichstatuses:
    useenrichstatuses = enrichstatuses
  def projections():
#    yield Projections("170203", "2e2mu", "fa3", "rescalemixtures", "fullrange", "VHHadrtagged")
#    return
#    for production in config.productionsforcombine:
      production = None
      for analysis in analyses:
        if analysis not in useanalyses: continue
        for normalization in normalizations:
          for enrichstatus in enrichstatuses:
            if enrichstatus == "blind": continue
            if enrichstatus not in useenrichstatuses: continue
            yield Projections(analysis, normalization, production, enrichstatus)

  length = len(list(projections()))
  for i, (p, ch, ca) in enumerate(itertools.product(projections(), channels, categories), start=1):
    process = int(sys.argv[1])
    if process == 1 or process == 4:
      p.projections(ch, ca, nicestyle=True)
#      p.projections(ch, ca, nicestyle=True, Dbkg_allcategories=True)
      p.projections(ch, ca, nicestyle=True, forWIN=True)
#      p.projections(ch, ca, nicestyle=True, Dbkg_allcategories=True, forWIN=True)
    if process == 2 or process == 5:
      pass
    if process == 3:
      p.projections(ch, ca)
      p.projections(ch, ca, subdir="ggH", productionmode="ggH")
      p.projections(ch, ca, subdir="VBF", productionmode="VBF")
      p.projections(ch, ca, subdir="VH",  productionmode="VH")
    if process == 4:
      p.animation(ca, ch, nicestyle=True)
      p.animation(ca, ch, nicestyle=True, Dbkg_allcategories=True)
      p.animation(ca, ch, nicestyle=True, forWIN=True)
      p.animation(ca, ch, nicestyle=True, Dbkg_allcategories=True, forWIN=True)
    if process == 5:
      pass
    if process == 6:
      if p.enrichstatus == "fullrange":
        p.animation(ca, ch)
    print i, "/", length*len(channels)*len(categories)

if __name__ == "__main__":
    main()
