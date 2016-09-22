import abc
import collections
from helperstuff import config, constants
from helperstuff.combinehelpers import getrate, gettemplate
from helperstuff.enums import analyses, Analysis, categories, Category, Channel, channels, EnumItem, MultiEnum, MultiEnumABCMeta, MyEnum, Production, productions, Systematic, WhichProdDiscriminants
from helperstuff.filemanager import tfiles
from helperstuff.samples import ReweightingSample
import helperstuff.style
from helperstuff.templates import IntTemplate, Template, TemplatesFile
import itertools
from math import sqrt
import os
import ROOT
import subprocess

def GetAxis(self, i):
    if i == 0:
        return self.GetXaxis()
    elif i == 1:
        return self.GetYaxis()
    elif i == 2:
        return self.GetZaxis()
    assert False
ROOT.TH1.GetAxis = GetAxis
del GetAxis

def Projection(self, i):
    if i == 0:
        return self.ProjectionX()
    elif i == 1:
        return self.ProjectionY()
    elif i == 2:
        return self.ProjectionZ()
    assert False
ROOT.TH1.Projection = Projection
del Projection

class TemplateForProjection(object):
    __metaclass__ = abc.ABCMeta
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

class Normalization(MyEnum):
    enumname = "normalization"
    enumitems = [
                 EnumItem("defaultnormalization"),
                 EnumItem("rescalemixtures"),
                 EnumItem("areanormalize"),
                ]
normalizations = Normalization.items()

class EnrichStatus(MyEnum):
    enumname = "enrichstatus"
    enumitems = [
                 EnumItem("enrich"),
                 EnumItem("noenrich"),
                 EnumItem("impoverish"),
                ]
    def cuttext(self):
        if self == "enrich": return "D_{bkg} > 0.5"
        if self == "noenrich": return ""
        if self == "impoverish": return "D_{bkg} < 0.5"
        assert False
    def dirname(self):
        if self == "enrich": return "enrich"
        if self == "noenrich": return "fullrange"
        if self == "impoverish": return "blind"
        assert False
enrichstatuses = EnrichStatus.items()

class BaseTemplateFromFile(TemplateForProjection):
    def __init__(self, color, *args):
        super(BaseTemplateFromFile, self).__init__(*args)
        self.color = color
        self.h = self.template.gettemplate().Clone()
        if self.productionmode in ["ggH", "VBF"]:
            scalefactor = getrate("2e2mu", self.category, self.productionmode, "fordata", self.production) / Template(self.production, self.category, self.analysis, "2e2mu", self.productionmode, "0+", "prod+dec", self.whichproddiscriminants).gettemplate().Integral()
        elif self.productionmode == "data":
            scalefactor = 1
        else:
            scalefactor = getrate(self.channel, self.category, self.productionmode, "fordata", self.production) / self.Integral()
        for x, y, z in itertools.product(range(1, self.h.GetNbinsX()+1), range(1, self.h.GetNbinsY()+1), range(1, self.h.GetNbinsZ()+1)):
            if (self.enrichstatus == "impoverish" and self.h.GetZaxis().GetBinLowEdge(z) >= .5 \
             or self.enrichstatus == "enrich"     and self.h.GetZaxis().GetBinLowEdge(z) < .5):
                self.h.SetBinContent(x, y, z, 0)
        self.projections = {}
        if self.productionmode == "data":
            self.hstackoption = "P"
            self.h.SetMarkerColor(self.color)
            self.h.SetMarkerStyle(20)
        else:
            self.hstackoption = "hist"
        self.Scale(scalefactor)
    def check(self, *args):
        if self.normalization is None:
            self.normalization = Normalization("defaultnormalization")
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
    def __init__(self, title, color, SMintegral, *templatesandfactors):
        self.projections = {}
        self.title = title
        self.color = color
        self.h = None
        self.templatesandfactors = templatesandfactors
        for template, factor in templatesandfactors:
            if self.h is None:
                self.h = template.h.Clone()
                self.h.Scale(factor)
            else:
                self.h.Add(template.h, factor)
        analyses = set(t[0].analysis for t in templatesandfactors)
        assert len(analyses) == 1
        analyses = list(analyses)[0]
        if normalization == "rescalemixtures":
            self.Scale(SMintegral / self.Integral())
        else:
            self.Scale(.5)
        self.hstackoption = "hist"

    @property
    def discriminants(self):
        assert len({tf[0].discriminants for tf in self.templatesandfactors}) == 1
        return self.templatesandfactors[0][0].discriminants

exts = "png", "eps", "root", "pdf"


class Projections(MultiEnum):
  enums = [Channel, Analysis, Normalization, Systematic, Production, EnrichStatus, Category, WhichProdDiscriminants]
  enumname = "projections"

  def applysynonyms(self, enumsdict):
    if enumsdict[WhichProdDiscriminants] is None:
      enumsdict[WhichProdDiscriminants] = "D_g13gi1"
    super(Projections, self).applysynonyms(enumsdict)

  def check(self, *args):
    if self.normalization is None:
      self.normalization = Normalization("defaultnormalization")
    if self.systematic is None:
      self.systematic = Systematic("")
    super(Projections, self).check(*args)

  def projections(self):
    print self

    giname = self.analysis.couplingname

    ggHSM     = TemplateFromFile(   1,              "ggH", "prod+dec", self.category, self.whichproddiscriminants, self.enrichstatus, self.normalization, self.production, self.channel, self.systematic, self.analysis, self.analysis.purehypotheses[0])
    ggHBSM    = TemplateFromFile(   ROOT.kCyan,     "ggH", "prod+dec", self.category, self.whichproddiscriminants, self.enrichstatus, self.normalization, self.production, self.channel, self.systematic, self.analysis, self.analysis.purehypotheses[1])
    ggHint    = IntTemplateFromFile(0,              "ggH", "prod+dec", self.category, self.whichproddiscriminants, self.enrichstatus, self.normalization, self.production, self.channel, self.systematic, self.analysis, "g11gi1")

    ggHmix_p  = TemplateSum("ggH {}=#plus0.5" .format(self.analysis.title), ROOT.kGreen+3, ggHSM.Integral(), (ggHSM, 1), (ggHBSM, 1), (ggHint,  1))
    ggHmix_m  = TemplateSum("ggH {}=#minus0.5".format(self.analysis.title), 4,             ggHSM.Integral(), (ggHSM, 1), (ggHBSM, 1), (ggHint, -1))

    VBFSM = \
    VBFg14gi0 = TemplateFromFile(   ROOT.kAzure+1,  "VBF", "prod+dec", self.category, self.whichproddiscriminants, self.enrichstatus, self.normalization, self.production, self.channel, self.systematic, self.analysis, self.analysis.purehypotheses[0])
    VBFg13gi1 = IntTemplateFromFile(0,              "VBF", "prod+dec", self.category, self.whichproddiscriminants, self.enrichstatus, self.normalization, self.production, self.channel, self.systematic, self.analysis, "g13gi1")
    VBFg12gi2 = IntTemplateFromFile(0,              "VBF", "prod+dec", self.category, self.whichproddiscriminants, self.enrichstatus, self.normalization, self.production, self.channel, self.systematic, self.analysis, "g12gi2")
    VBFg11gi3 = IntTemplateFromFile(0,              "VBF", "prod+dec", self.category, self.whichproddiscriminants, self.enrichstatus, self.normalization, self.production, self.channel, self.systematic, self.analysis, "g11gi3")
    VBFg10gi4 = TemplateFromFile(   0,              "VBF", "prod+dec", self.category, self.whichproddiscriminants, self.enrichstatus, self.normalization, self.production, self.channel, self.systematic, self.analysis, self.analysis.purehypotheses[1])

    VBFpieces = [VBFg14gi0, VBFg13gi1, VBFg12gi2, VBFg11gi3, VBFg10gi4]

    gi_VBFBSM = (ReweightingSample("VBF", "SM").xsec / ReweightingSample("VBF", self.analysis.purehypotheses[1]).xsec)**.25
    g1_VBFmix = sqrt(2)
    gi_VBFmix = sqrt(2)*gi_BSM

    VBFBSM    = TemplateSum(VBFg10gi4.title,                                3,             VBFSM.Integral(), (VBFg10gi4, gi_BSM**4))

    VBFmix_p  = TemplateSum("VBF {}=#plus0.5" .format(self.analysis.title), 5,            VBFSM.Integral(),
                            *((template, g1_mix**(4-j) * gi_mix**j) for j, template in enumerate(VBFpieces))
                           )
    VBFmix_m  = TemplateSum("VBF {}=#minus0.5".format(self.analysis.title), ROOT.kTeal-1,  VBFSM.Integral(),
                            *((template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(VBFpieces))
                           )

    ZHSM = \
    ZHg14gi0 = TemplateFromFile(   ROOT.kAzure+1,  "ZH", "prod+dec", self.category, self.whichproddiscriminants, self.enrichstatus, self.normalization, self.production, self.channel, self.systematic, self.analysis, self.analysis.purehypotheses[0])
    ZHg13gi1 = IntTemplateFromFile(0,              "ZH", "prod+dec", self.category, self.whichproddiscriminants, self.enrichstatus, self.normalization, self.production, self.channel, self.systematic, self.analysis, "g13gi1")
    ZHg12gi2 = IntTemplateFromFile(0,              "ZH", "prod+dec", self.category, self.whichproddiscriminants, self.enrichstatus, self.normalization, self.production, self.channel, self.systematic, self.analysis, "g12gi2")
    ZHg11gi3 = IntTemplateFromFile(0,              "ZH", "prod+dec", self.category, self.whichproddiscriminants, self.enrichstatus, self.normalization, self.production, self.channel, self.systematic, self.analysis, "g11gi3")
    ZHg10gi4 = TemplateFromFile(   0,              "ZH", "prod+dec", self.category, self.whichproddiscriminants, self.enrichstatus, self.normalization, self.production, self.channel, self.systematic, self.analysis, self.analysis.purehypotheses[1])

    ZHpieces = [ZHg14gi0, ZHg13gi1, ZHg12gi2, ZHg11gi3, ZHg10gi4]

    gi_ZHBSM = (ReweightingSample("ZH", "SM").xsec / ReweightingSample("ZH", self.analysis.purehypotheses[1]).xsec)**.25
    g1_ZHmix = sqrt(2)
    gi_ZHmix = sqrt(2)*gi_BSM

    ZHBSM    = TemplateSum(ZHg10gi4.title,                                3,             ZHSM.Integral(), (ZHg10gi4, gi_BSM**4))

    ZHmix_p  = TemplateSum("ZH {}=#plus0.5" .format(self.analysis.title), 5,            ZHSM.Integral(),
                            *((template, g1_mix**(4-j) * gi_mix**j) for j, template in enumerate(ZHpieces))
                           )
    ZHmix_m  = TemplateSum("ZH {}=#minus0.5".format(self.analysis.title), ROOT.kTeal-1,  ZHSM.Integral(),
                            *((template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(ZHpieces))
                           )

    WHSM = \
    WHg14gi0 = TemplateFromFile(   ROOT.kAzure+1,  "WH", "prod+dec", self.category, self.whichproddiscriminants, self.enrichstatus, self.normalization, self.production, self.channel, self.systematic, self.analysis, self.analysis.purehypotheses[0])
    WHg13gi1 = IntTemplateFromFile(0,              "WH", "prod+dec", self.category, self.whichproddiscriminants, self.enrichstatus, self.normalization, self.production, self.channel, self.systematic, self.analysis, "g13gi1")
    WHg12gi2 = IntTemplateFromFile(0,              "WH", "prod+dec", self.category, self.whichproddiscriminants, self.enrichstatus, self.normalization, self.production, self.channel, self.systematic, self.analysis, "g12gi2")
    WHg11gi3 = IntTemplateFromFile(0,              "WH", "prod+dec", self.category, self.whichproddiscriminants, self.enrichstatus, self.normalization, self.production, self.channel, self.systematic, self.analysis, "g11gi3")
    WHg10gi4 = TemplateFromFile(   0,              "WH", "prod+dec", self.category, self.whichproddiscriminants, self.enrichstatus, self.normalization, self.production, self.channel, self.systematic, self.analysis, self.analysis.purehypotheses[1])

    WHpieces = [WHg14gi0, WHg13gi1, WHg12gi2, WHg11gi3, WHg10gi4]

    gi_WHBSM = (ReweightingSample("WH", "SM").xsec / ReweightingSample("WH", self.analysis.purehypotheses[1]).xsec)**.25
    g1_WHmix = sqrt(2)
    gi_WHmix = sqrt(2)*gi_BSM

    WHBSM    = TemplateSum(WHg10gi4.title,                                3,             WHSM.Integral(), (WHg10gi4, gi_BSM**4))

    WHmix_p  = TemplateSum("WH {}=#plus0.5" .format(self.analysis.title), 5,            WHSM.Integral(),
                            *((template, g1_mix**(4-j) * gi_mix**j) for j, template in enumerate(WHpieces))
                           )
    WHmix_m  = TemplateSum("WH {}=#minus0.5".format(self.analysis.title), ROOT.kTeal-1,  WHSM.Integral(),
                            *((template, g1_mix**(4-j) * (-gi_mix)**j) for j, template in enumerate(WHpieces))
                           )

    qqZZ      = TemplateFromFile(6,              "prod+dec", self.category, self.whichproddiscriminants, self.enrichstatus, self.normalization, self.analysis, self.channel, "qqZZ",    self.systematic, self.production)
    ggZZ      = TemplateFromFile(ROOT.kOrange+6, "prod+dec", self.category, self.whichproddiscriminants, self.enrichstatus, self.normalization, self.analysis, self.channel, "ggZZ",    self.systematic, self.production)
    VBFbkg    = TemplateFromFile(ROOT.kViolet+7, "prod+dec", self.category, self.whichproddiscriminants, self.enrichstatus, self.normalization, self.analysis, self.channel, "VBF bkg", self.systematic, self.production)
    ZX        = TemplateFromFile(2,              "prod+dec", self.category, self.whichproddiscriminants, self.enrichstatus, self.normalization, self.analysis, self.channel, "ZX",      self.systematic, self.production)

    templates = [
                 ggHSM, ggHBSM, ggHint, ggHmix_p, ggHmix_m,
                ] + VBFpieces + [
                 VBFBSM, VBFmix_p, VBFmix_m,
                ] + ZHpieces + [
                 ZHBSM, ZHmix_p, ZHmix_m,
                ] + WHpieces + [
                 WHBSM, WHmix_p, WHmix_m,
                ] + [
                 qqZZ, ggZZ, VBFbkg, ZX
                ]

    if self.enrichstatus == "impoverish" and config.usedata or config.unblinddistributions:
        templates += [
                      TemplateFromFile(1, "prod+dec", self.category, self.whichproddiscriminants, self.enrichstatus, self.normalization, self.analysis, self.channel, "data", self.production, "unblind" if config.unblinddistributions else "blind")
                     ]

    c1 = ROOT.TCanvas()
    legend = ROOT.TLegend(.75, .65, .9, .9)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    for template in templates:
        if template.color:
            if self.normalization == "areanormalize":
                try:
                    template.Scale(1/template.Integral())
                except ZeroDivisionError:
                    pass
            template.AddToLegend(legend)

    for i, discriminant in enumerate(TemplatesFile(self.channel, self.systematic, "ggh", self.analysis, self.production, "prod+dec", self.whichproddiscriminants, self.category).discriminants):
        hstack = ROOT.THStack(discriminant.name, discriminant.title)
        for template in templates:
            if template.color:
                hstack.Add(template.Projection(i), template.hstackoption)
        hstack.Draw("nostack")
        hstack.GetXaxis().SetTitle(discriminant.title)
        legend.Draw()
        try:
            os.makedirs(self.saveasdir)
        except OSError:
            pass
        for ext in exts:
            c1.SaveAs(os.path.join(self.saveasdir, "{}.{}".format(discriminant.name, ext)))

  @property
  def saveasdir(self):
      return os.path.join(config.plotsbasedir, "templateprojections", self.enrichstatus.dirname(), str(self.normalization), "{}_{}/{}/{}".format(self.analysis, self.production, self.category, self.channel))


def projections(*args):
    Projections(*args).projections()

if __name__ == "__main__":
  for production in productions:
    for channel in channels:
      for analysis in analyses:
        for normalization in normalizations:
          for enrichstatus in enrichstatuses:
            for category in categories:
              if normalization != "rescalemixtures": continue   #uncomment this to get the niceplots fast
#              if channel != "2e2mu" or analysis != "fa3" or normalization != "areanormalize": continue
              print production, channel, analysis, normalization, enrichstatus, category
              projections(channel, analysis, normalization, production, enrichstatus, category)
