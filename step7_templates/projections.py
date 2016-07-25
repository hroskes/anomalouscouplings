import collections
from helperstuff import config, constants
from helperstuff.combinehelpers import getrate
from helperstuff.enums import analyses, Analysis, Channel, channels, EnumItem, MultiEnum, MyEnum, Production, productions, Systematic, Template
from helperstuff.filemanager import tfiles
import helperstuff.style
import itertools
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

Axis = collections.namedtuple("Axis", "name title index")
axes = [
        None,
        None,
        Axis("Dbkg", "D_{bkg}", 2),
       ]

class TemplateForProjection(object):
    def Projection(self, i):
        if i not in self.projections:
            result = self.h.Projection(i)
            result.SetName(axes[i].name)
            result.SetTitle(axes[i].title)
            result.SetXTitle(axes[i].title)
            result.SetLineColor(self.color)
            self.projections[i] = result
        return self.projections[i]

    def AddToLegend(self, legend):
        legend.AddEntry(self.Projection(0), self.title, "l")

    def Scale(self, *args, **kwargs):
        return self.h.Scale(*args, **kwargs)

    def Integral(self, *args, **kwargs):
        return self.h.Integral(*args, **kwargs)

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
        if self == "enrich": return "D_{bkg}>0.5"
        if self == "noenrich": return ""
        if self == "impoverish": return "D_{bkg}<0.5"
        assert False
    def dirname(self):
        if self == "enrich": return "enrich"
        if self == "noenrich": return "fullrange"
        if self == "impoverish": return "blind"
        assert False
enrichstatuses = EnrichStatus.items()
    
class TemplateFromFile(TemplateForProjection, MultiEnum):
    enums = [Template, Normalization, EnrichStatus]
    def __init__(self, color, *args):
        super(TemplateFromFile, self).__init__(*args)
        self.title = self.template.title()
        self.color = color
        self.h = self.template.gettemplate().Clone()
        if self.productionmode == "ggH":
            scalefactor = getrate("2e2mu", self.productionmode, "fordata", self.production) / Template(self.production, self.analysis, "2e2mu", self.productionmode, "0+").gettemplate().Integral()
        elif self.productionmode == "data":
            scalefactor = 1
        else:
            scalefactor = getrate(self.channel, self.productionmode, "fordata", self.production) / self.Integral()
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
        super(TemplateFromFile, self).check(*args)

class TemplateSum(TemplateForProjection):
    def __init__(self, title, color, SMintegral, *templatesandfactors):
        self.projections = {}
        self.title = title
        self.color = color
        self.h = None
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

exts = "png", "eps", "root", "pdf"


class Projections(MultiEnum):
  enums = [Channel, Analysis, Normalization, Systematic, Production, EnrichStatus]
  enumname = "projections"
  def check(self, *args):
    if self.normalization is None:
      self.normalization = Normalization("defaultnormalization")
    if self.systematic is None:
      self.systematic = Systematic("")
    super(Projections, self).check(*args)

  def projections(self):
    print self
    templates = [
                 TemplateFromFile(1, self.enrichstatus, self.normalization, self.analysis.signaltemplates(self.production, self.channel, self.systematic)[0]),
                 TemplateFromFile(ROOT.kCyan, self.enrichstatus, self.normalization, self.analysis.signaltemplates(self.production, self.channel, self.systematic)[1]),
                 TemplateFromFile(0, self.enrichstatus, self.normalization, self.analysis.signaltemplates(self.production, self.channel, self.systematic)[2]),
                ]
    templates+= [
                 TemplateSum("ggH {}=0.5".format(self.analysis.title), ROOT.kGreen+3, templates[0].Integral(), (templates[0], 1), (templates[1], 1), (templates[2], 1)),
                 TemplateSum("ggH {}=-0.5".format(self.analysis.title), 4, templates[0].Integral(), (templates[0], 1), (templates[1], 1), (templates[2], -1)),
                 TemplateFromFile(6, self.enrichstatus, self.normalization, self.analysis, self.channel, "qqZZ", self.systematic, self.production),
                 TemplateFromFile(ROOT.kOrange+6, self.enrichstatus, self.normalization, self.analysis, self.channel, "ggZZ", self.systematic, self.production),
                 TemplateFromFile(2, self.enrichstatus, self.normalization, self.analysis, self.channel, "ZX", self.systematic, self.production),
                ]
    if self.enrichstatus == "impoverish" or config.unblinddistributions:
        templates += [
                      TemplateFromFile(1, self.enrichstatus, self.normalization, self.analysis, self.channel, "data", self.production, "unblind" if config.unblinddistributions else "blind")
                     ]
    axes[0] = Axis(self.analysis.purediscriminant(), self.analysis.purediscriminant(title=True), 0)
    axes[1] = Axis(self.analysis.mixdiscriminant(), self.analysis.mixdiscriminant(title=True), 1)

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

    for axis in axes:
        hstack = ROOT.THStack(axis.name, axis.title)
        for template in templates:
            if template.color:
                hstack.Add(template.Projection(axis.index), template.hstackoption)
        hstack.Draw("nostack")
        hstack.GetXaxis().SetTitle(axis.title)
        legend.Draw()
        try:
            os.makedirs(self.saveasdir)
        except OSError:
            pass
        for ext in exts:
            c1.SaveAs(os.path.join(self.saveasdir, "{}.{}".format(axis.name, ext)))

  @property
  def saveasdir(self):
      return os.path.join(config.plotsbasedir, "templateprojections", self.enrichstatus.dirname(), str(self.normalization), "{}_{}/{}".format(self.analysis, self.production, self.channel))


def projections(*args):
    Projections(*args).projections()

if __name__ == "__main__":
  for production in productions:
    for channel in channels:
      for analysis in analyses:
        for normalization in normalizations:
          for enrichstatus in enrichstatuses:
#            if channel != "2e2mu" or analysis != "fa3" or normalization != "areanormalize": continue
            projections(channel, analysis, normalization, production, enrichstatus)
