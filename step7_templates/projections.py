#!/usr/bin/env python
import abc
import collections
from helperstuff import config, constants, run1info
from helperstuff.combinehelpers import getrate, gettemplate
from helperstuff.enums import analyses, Analysis, categories, Category, Channel, channels, EnumItem, MultiEnum, MultiEnumABCMeta, MyEnum, Production, ProductionMode, productions, ShapeSystematic
from helperstuff.samples import ReweightingSample, samplewithfai
import helperstuff.rootoverloads.histogramaxisnumbers
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
        return self.h.Scale(*args, **kwargs)

    def Integral(self, *args, **kwargs):
        return self.h.Integral(*args, **kwargs)

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
        return "h", "hstackoption"

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
        super(BaseTemplateFromFile, self).__init__(*args)
        self.color = color
        #I don't know what the following line accomplishes
        # but it fixes the WH interference templates.
        #Otherwise they are wrong in some way that
        # I don't understand
        self.template.gettemplate()

    def inithistogram(self):
        self.h = self.template.gettemplate().Clone()

        if self.productionmode in ["ggH", "VBF", "ZH", "WH"]:
            scalefactor = getrate("2e2mu", self.category, self.productionmode, "fordata", self.production, self.analysis) / Template(self.production, self.category, self.analysis, "2e2mu", self.productionmode, self.analysis.purehypotheses[0]).gettemplate().Integral()
        elif self.productionmode == "data" or self.Integral() == 0:
            scalefactor = 1
        else:
            scalefactor = getrate(self.channel, self.category, self.productionmode, "fordata", self.production, self.analysis) / self.Integral()

        for x, y, z in itertools.product(range(1, self.h.GetNbinsX()+1), range(1, self.h.GetNbinsY()+1), range(1, self.h.GetNbinsZ()+1)):
            if (self.enrichstatus == "impoverish" and self.h.GetZaxis().GetBinLowEdge(z) >= .5 \
             or self.enrichstatus == "enrich"     and self.h.GetZaxis().GetBinLowEdge(z) < .5):
                self.h.SetBinContent(x, y, z, 0)
        if self.productionmode == "data":
            self.hstackoption = "P"
            self.h.SetMarkerColor(self.color)
            self.h.SetMarkerStyle(20)
        else:
            self.hstackoption = "hist"

        self.Scale(scalefactor)

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

class TemplateSum(TemplateForProjection):
    def __init__(self, title, color, *templatesandfactors):
        self.title = title
        self.color = color
        self.templatesandfactors = templatesandfactors
        analyses = {t[0].analysis for t in templatesandfactors}
        assert len(analyses) == 1
        self.analysis = analyses.pop()
        super(TemplateSum, self).__init__()

    def inithistogram(self):
        self.h = None
        for template, factor in self.templatesandfactors:
            if not factor: continue
            if self.h is None:
                self.h = template.h.Clone()
                self.h.Scale(factor)
            else:
                self.h.Add(template.h, factor)
        self.hstackoption = "hist"
        super(TemplateSum, self).inithistogram()

    @property
    def discriminants(self):
        assert len({template.discriminants for template, factor in self.templatesandfactors}) == 1
        return self.templatesandfactors[0][0].discriminants

class ComponentTemplateSum(TemplateSum):
    def __init__(self, title, color, SMintegral, *templatesandfactors):
        """
        Works the same as TemplateSum, but rescales
        """
        self.SMintegral = SMintegral
        super(ComponentTemplateSum, self).__init__(title, color, *templatesandfactors)
        for template, factor in templatesandfactors:
            if not isinstance(template, BaseTemplateFromFile):
                raise TypeError("ComponentTemplateSum can only come from TemplatesFromFiles")

    def inithistogram(self):
        super(ComponentTemplateSum, self).inithistogram()
        if self.Integral():
            self.Scale(self.SMintegral / self.Integral())

class Projections(MultiEnum):
  enums = [Channel, Analysis, Normalization, ShapeSystematic, Production, EnrichStatus, Category]
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

  def projections(self, **kwargs):
    giname = self.analysis.couplingname
    BSMhypothesis = self.analysis.purehypotheses[1]

    ggHfactor = VBFfactor = VHfactor = 1
    subdir = ""
    saveasdir = self.saveasdir
    customfai = None
    animation = False
    saveasappend = ""
    exts = "png", "eps", "root", "pdf"
    otherthingstodraw = []
    for kw, kwarg in kwargs.iteritems():
       if kw == "ggHfactor":
           ggHfactor = float(kwarg)
       elif kw == "VBFfactor":
           VBFfactor = float(kwarg)
       elif kw == "VHfactor":
           VHfactor = float(kwarg)
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
       else:
           raise TypeError("Unknown kwarg {}={}".format(kw, kwarg))

    SMhypothesis = self.analysis.purehypotheses[0]
    gi_ggHBSM = getattr(ReweightingSample("ggH", BSMhypothesis), BSMhypothesis.couplingname)
    gi_VBFBSM = copysign((ReweightingSample("VBF", SMhypothesis).xsec / ReweightingSample("VBF", BSMhypothesis).xsec)**.25, gi_ggHBSM)
    gi_VHBSM = copysign(((ReweightingSample("WH", SMhypothesis).xsec + ReweightingSample("ZH", SMhypothesis).xsec) / (ReweightingSample("WH", BSMhypothesis).xsec + ReweightingSample("ZH", BSMhypothesis).xsec))**.25, gi_ggHBSM)
    if self.category == "UntaggedIchep16":
        g1_mix = 1/sqrt(2)
        gi_mix = 1/sqrt(2)*gi_ggHBSM
        fainame = "{}^{{{}}}".format(self.analysis.title, "dec")
    elif self.category == "VBF2jTaggedIchep16":
        g1_mix = 1/2**.25
        gi_mix = 1/2**.25 * gi_VBFBSM
        fainame = "{}^{{{}}}".format(self.analysis.title, "VBFdec")
    elif self.category == "VHHadrTaggedIchep16":
        g1_mix = 1/2**.25
        gi_mix = 1/2**.25 * gi_VHBSM
        fainame = "{}^{{{}}}".format(self.analysis.title, "VHdec")

    ggHSM     = self.TemplateFromFile(   0, "ggH", self.category, self.enrichstatus, self.normalization, self.production, self.channel, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
    ggHBSM    = self.TemplateFromFile(   0, "ggH", self.category, self.enrichstatus, self.normalization, self.production, self.channel, self.shapesystematic, self.analysis, BSMhypothesis)
    ggHint    = self.IntTemplateFromFile(0, "ggH", self.category, self.enrichstatus, self.normalization, self.production, self.channel, self.shapesystematic, self.analysis, "g11gi1")

    ggHmix_p  = ComponentTemplateSum("ggH {}=#plus0.5" .format(fainame), 0, ggHSM.Integral(), (ggHSM, g1_mix**2), (ggHBSM, (gi_mix/gi_ggHBSM)**2), (ggHint,  g1_mix*gi_mix/gi_ggHBSM))
    ggHmix_m  = ComponentTemplateSum("ggH {}=#minus0.5".format(fainame), 0, ggHSM.Integral(), (ggHSM, g1_mix**2), (ggHBSM, (gi_mix/gi_ggHBSM)**2), (ggHint,  -g1_mix*gi_mix/gi_ggHBSM))

    VBFSM = \
    VBFg14gi0 = self.TemplateFromFile(   0, "VBF", self.category, self.enrichstatus, self.normalization, self.production, self.channel, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
    VBFg13gi1 = self.IntTemplateFromFile(0, "VBF", self.category, self.enrichstatus, self.normalization, self.production, self.channel, self.shapesystematic, self.analysis, "g13gi1")
    VBFg12gi2 = self.IntTemplateFromFile(0, "VBF", self.category, self.enrichstatus, self.normalization, self.production, self.channel, self.shapesystematic, self.analysis, "g12gi2")
    VBFg11gi3 = self.IntTemplateFromFile(0, "VBF", self.category, self.enrichstatus, self.normalization, self.production, self.channel, self.shapesystematic, self.analysis, "g11gi3")
    VBFg10gi4 = self.TemplateFromFile(   0, "VBF", self.category, self.enrichstatus, self.normalization, self.production, self.channel, self.shapesystematic, self.analysis, BSMhypothesis)

    VBFpieces = [VBFg14gi0, VBFg13gi1, VBFg12gi2, VBFg11gi3, VBFg10gi4]

    VBFBSM    = ComponentTemplateSum(VBFg10gi4.title,                                0, VBFSM.Integral(), (VBFg10gi4, gi_VBFBSM**4))

    VBFmix_p  = ComponentTemplateSum("VBF {}=#plus0.5" .format(fainame), 0, VBFSM.Integral(),
                                     *((template, g1_mix**(4-j) * gi_mix**j) for j, template in enumerate(VBFpieces))
                                    )
    VBFmix_m  = ComponentTemplateSum("VBF {}=#minus0.5".format(fainame), 0, VBFSM.Integral(),
                                     *((template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(VBFpieces))
                                    )

    ZHSM = \
    ZHg14gi0 = self.TemplateFromFile(   0, "ZH", self.category, self.enrichstatus, self.normalization, self.production, self.channel, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
    ZHg13gi1 = self.IntTemplateFromFile(0, "ZH", self.category, self.enrichstatus, self.normalization, self.production, self.channel, self.shapesystematic, self.analysis, "g13gi1")
    ZHg12gi2 = self.IntTemplateFromFile(0, "ZH", self.category, self.enrichstatus, self.normalization, self.production, self.channel, self.shapesystematic, self.analysis, "g12gi2")
    ZHg11gi3 = self.IntTemplateFromFile(0, "ZH", self.category, self.enrichstatus, self.normalization, self.production, self.channel, self.shapesystematic, self.analysis, "g11gi3")
    ZHg10gi4 = self.TemplateFromFile(   0, "ZH", self.category, self.enrichstatus, self.normalization, self.production, self.channel, self.shapesystematic, self.analysis, BSMhypothesis)

    ZHpieces = [ZHg14gi0, ZHg13gi1, ZHg12gi2, ZHg11gi3, ZHg10gi4]

    ZHBSM    = ComponentTemplateSum(ZHg10gi4.title,                                0, ZHSM.Integral(), (ZHg10gi4, gi_VHBSM**4))

    ZHmix_p  = ComponentTemplateSum("ZH {}=#plus0.5" .format(fainame), 0, ZHSM.Integral(),
                                    *((template, g1_mix**(4-j) * gi_mix**j) for j, template in enumerate(ZHpieces))
                                   )
    ZHmix_m  = ComponentTemplateSum("ZH {}=#minus0.5".format(fainame), 0, ZHSM.Integral(),
                                    *((template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(ZHpieces))
                                   )

    WHSM = \
    WHg14gi0 = self.TemplateFromFile(   0, "WH", self.category, self.enrichstatus, self.normalization, self.production, self.channel, self.shapesystematic, self.analysis, self.analysis.purehypotheses[0])
    WHg13gi1 = self.IntTemplateFromFile(0, "WH", self.category, self.enrichstatus, self.normalization, self.production, self.channel, self.shapesystematic, self.analysis, "g13gi1")
    WHg12gi2 = self.IntTemplateFromFile(0, "WH", self.category, self.enrichstatus, self.normalization, self.production, self.channel, self.shapesystematic, self.analysis, "g12gi2")
    WHg11gi3 = self.IntTemplateFromFile(0, "WH", self.category, self.enrichstatus, self.normalization, self.production, self.channel, self.shapesystematic, self.analysis, "g11gi3")
    WHg10gi4 = self.TemplateFromFile(   0, "WH", self.category, self.enrichstatus, self.normalization, self.production, self.channel, self.shapesystematic, self.analysis, BSMhypothesis)

    WHpieces = [WHg14gi0, WHg13gi1, WHg12gi2, WHg11gi3, WHg10gi4]

    WHBSM    = ComponentTemplateSum(WHg10gi4.title,                                0, WHSM.Integral(), (WHg10gi4, gi_VHBSM**4))

    WHmix_p  = ComponentTemplateSum("WH {}=#plus0.5" .format(fainame), 0, WHSM.Integral(),
                                    *((template, g1_mix**(4-j) * gi_mix**j) for j, template in enumerate(WHpieces))
                                   )
    WHmix_m  = ComponentTemplateSum("WH {}=#minus0.5".format(fainame), 0, WHSM.Integral(),
                                    *((template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(WHpieces))
                                   )

    SM    = TemplateSum("SM",                                1,             (ggHSM,    ggHfactor), (VBFSM,    VBFfactor), (ZHSM,    VHfactor), (WHSM,    VHfactor))
    BSM   = TemplateSum("{}=1".format(self.analysis.title),  ROOT.kCyan,    (ggHBSM,   ggHfactor), (VBFBSM,   VBFfactor), (ZHBSM,   VHfactor), (WHBSM,   VHfactor))
    mix_p = TemplateSum("{}=#plus0.5".format(fainame),       ROOT.kGreen+3, (ggHmix_p, ggHfactor), (VBFmix_p, VBFfactor), (ZHmix_p, VHfactor), (WHmix_p, VHfactor))
    mix_m = TemplateSum("{}=#minus0.5".format(fainame),      4,             (ggHmix_m, ggHfactor), (VBFmix_m, VBFfactor), (ZHmix_m, VHfactor), (WHmix_m, VHfactor))

    qqZZ      = self.TemplateFromFile(6,              self.category, self.enrichstatus, self.normalization, self.analysis, self.channel, "qqZZ",    self.shapesystematic, self.production)
    ggZZ      = self.TemplateFromFile(ROOT.kOrange+6, self.category, self.enrichstatus, self.normalization, self.analysis, self.channel, "ggZZ",    self.shapesystematic, self.production)
    VBFbkg    = self.TemplateFromFile(ROOT.kViolet+7, self.category, self.enrichstatus, self.normalization, self.analysis, self.channel, "VBF bkg", self.shapesystematic, self.production)
    ZX        = self.TemplateFromFile(2,              self.category, self.enrichstatus, self.normalization, self.analysis, self.channel, "ZX",      self.shapesystematic, self.production)

    templates = []
    if ggHfactor:
        templates += [
                      ggHSM, ggHBSM, ggHint, ggHmix_p, ggHmix_m,
                     ]
    if VBFfactor:
        templates += VBFpieces
        templates += [
                      VBFBSM, VBFmix_p, VBFmix_m,
                     ]
    if VHfactor:
        templates += ZHpieces
        templates += [
                      ZHBSM, ZHmix_p, ZHmix_m,
                     ]
        templates += WHpieces
        templates += [
                      WHBSM, WHmix_p, WHmix_m,
                     ]
    if not animation:
        templates += [
                      SM, BSM, mix_p, mix_m,
                     ]
    else:
        assert ggHfactor == VBFfactor == VHfactor == 1
        if customfai <  0: plusminus = "#minus"
        if customfai == 0: plusminus = ""
        if customfai  > 0: plusminus = "#plus"
        ggHcustom = ComponentTemplateSum("ggH ({}^{{{}}}={}{:.2f})".format(self.analysis.title, "dec", plusminus, abs(customsample.fai("ggH", self.analysis))), 1, ggHSM.Integral(), (ggHSM, g1_custom**2), (ggHBSM, (gi_custom/gi_ggHBSM)**2), (ggHint,  g1_custom*gi_custom/gi_ggHBSM))
        VBFcustom = ComponentTemplateSum("VBF ({}^{{{}}}={}{:.2f})".format(self.analysis.title, "VBF", plusminus, abs(customsample.fai("VBF", self.analysis))), 2, VBFSM.Integral(),
                                         *((template, g1_custom**(4-j) * gi_custom**j) for j, template in enumerate(VBFpieces))
                                        )
        ZHcustom  = ComponentTemplateSum("", 0, ZHSM.Integral(),
                                         *((template, g1_custom**(4-j) * gi_custom**j) for j, template in enumerate(ZHpieces))
                                        )
        WHcustom  = ComponentTemplateSum("", 0, WHSM.Integral(),
                                         *((template, g1_custom**(4-j) * gi_custom**j) for j, template in enumerate(WHpieces))
                                        )
        VHcustom = TemplateSum("VH ({}^{{{}}}={}{:.2f})".format(self.analysis.title, "VH", plusminus, abs(customsample.fai("VH", self.analysis))), 4, (ZHcustom, 1), (WHcustom, 1))
        for t in ggHcustom, VBFcustom, VHcustom:
            t.Scale(1.0/t.Integral())

        templates += [
                      ggHcustom, VBFcustom, ZHcustom, WHcustom, VHcustom
                     ]

    if ggHfactor == VBFfactor == VHfactor == 1 and not animation:
        templates += [
                      qqZZ, ggZZ, VBFbkg, ZX,
                     ]

    if self.enrichstatus == "impoverish" and config.usedata or config.unblinddistributions:
        templates += [
                      TemplateFromFile(1, self.category, self.enrichstatus, self.normalization, self.analysis, self.channel, "data", self.production, "unblind" if config.unblinddistributions else "blind")
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
            if template.title:
                template.AddToLegend(legend)

    for i, discriminant in enumerate(self.discriminants):
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

  @property
  def saveasdir(self):
      assert self.normalization == "rescalemixtures"
      return os.path.join(config.plotsbasedir, "templateprojections", self.enrichstatus.dirname(), "{}_{}/{}/{}".format(self.analysis, self.production, self.category, self.channel))

  @property
  def discriminants(self):
      return TemplatesFile(self.channel, self.shapesystematic, "ggh", self.analysis, self.production, self.category).discriminants

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

  def animation(self):
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
                    }

      for i, step in enumerate(animation):
          kwargs = kwargs_base.copy()
          kwargs["customfaiforanimation"] = step.fai_decay, "ggH"
          kwargs["saveasappend"] = i
          kwargs["otherthingstodraw"] = [step.excludedtext]
          self.projections(**kwargs)

      for discriminant in self.discriminants:
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
    #yield Projections("170119", "2e2mu", "fa3", "rescalemixtures", "fullrange", "Untagged")
    #yield Projections("170119", "2e2mu", "fa3", "rescalemixtures", "fullrange", "VBFtagged")
    #yield Projections("170119", "2e2mu", "fa3", "rescalemixtures", "fullrange", "VHHadrtagged")
    #return
    for production in productions:
      for channel in channels:
        for analysis in analyses:
          for normalization in normalizations:
            for enrichstatus in enrichstatuses:
              for category in categories:
                if normalization != "rescalemixtures": continue   #uncomment this to get the niceplots fast
                yield Projections(channel, analysis, normalization, production, enrichstatus, category)

  length = len(list(projections()))
  for i, p in enumerate(projections(), start=1):
    p.projections()
    p.projections(subdir="ggH", ggHfactor=1, VBFfactor=0, VHfactor=0)
    p.projections(subdir="VBF", ggHfactor=0, VBFfactor=1, VHfactor=0)
    p.projections(subdir="VH",  ggHfactor=0, VBFfactor=0, VHfactor=1)
    if p.enrichstatus == "fullrange":
      p.animation()
    print i, "/", length
