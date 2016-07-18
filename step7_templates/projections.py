import collections
from helperstuff import config, constants
from helperstuff.combinehelpers import getrate
from helperstuff.enums import analyses, Analysis, BlindStatus, blindstatuses, Channel, channels, EnumItem, MultiEnum, MyEnum, Production, productions, Systematic, Template
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
normalizations = [Normalization(normalization) for normalization in Normalization.enumitems]
    
class TemplateFromFile(TemplateForProjection, MultiEnum):
    enums = [Template, Normalization]
    def __init__(self, color, blindstatus, *args):  #have to do blindstatus separately, because Template has its own blindstatus :(
        super(TemplateFromFile, self).__init__(*args)
        self.blindstatus = BlindStatus(blindstatus)
        self.title = self.template.title()
        self.color = color
        self.h = self.template.gettemplate().Clone()
        if self.blindstatus == "blind":
            for x, y, z in itertools.product(range(1, self.h.GetNbinsX()+1), range(1, self.h.GetNbinsY()+1), range(1, self.h.GetNbinsZ()+1)):
                if self.h.GetZaxis().GetBinLowEdge(z) >= .5:
                    self.h.SetBinContent(x, y, z, 0)
        self.projections = {}
        if self.productionmode == "ggH":
            scalefactor = getrate("2e2mu", self.productionmode) / Template(self.production, self.analysis, "2e2mu", self.productionmode, "0+").gettemplate().Integral()
        elif self.productionmode == "data":
            scalefactor = 1
        else:
            scalefactor = getrate(self.channel, self.productionmode) / self.Integral()
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
    def __init__(self, title, color, mixturesign, *templatesandfactors):
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
        assert abs(mixturesign) == 1
        if normalization != "rescalemixtures":
            mixturesign = 0
        self.Scale(constants.JHUXS2L2la1 / (2*constants.JHUXS2L2la1 + mixturesign*analysis.interfxsec()))
        self.hstackoption = "hist"

exts = "png", "eps", "root", "pdf"


class Projections(MultiEnum):
  enums = [Channel, Analysis, Normalization, Systematic, Production, BlindStatus]
  enumname = "projections"
  def check(self, *args):
    if self.normalization is None:
      self.normalization = Normalization("defaultnormalization")
    if self.systematic is None:
      self.systematic = Systematic("")
    super(Projections, self).check(*args)

  def projections(self):
    templates = [
                 TemplateFromFile(1, self.blindstatus, self.normalization, self.analysis.signaltemplates(self.production, self.channel, self.systematic)[0]),
                 TemplateFromFile(ROOT.kCyan, self.blindstatus, self.normalization, self.analysis.signaltemplates(self.production, self.channel, self.systematic)[1]),
                 TemplateFromFile(0, self.blindstatus, self.normalization, self.analysis.signaltemplates(self.production, self.channel, self.systematic)[2]),
                ]
    templates+= [
                 TemplateSum("ggH {}=0.5".format(self.analysis.title()), ROOT.kGreen+3, 1, (templates[0], 1), (templates[1], 1), (templates[2], 1)),
                 TemplateSum("ggH {}=-0.5".format(self.analysis.title()), 4, -1, (templates[0], 1), (templates[1], 1), (templates[2], -1)),
                 TemplateFromFile(6, self.blindstatus, self.normalization, self.analysis, self.channel, "qqZZ", self.systematic, self.production),
                 TemplateFromFile(ROOT.kOrange+6, self.blindstatus, self.normalization, self.analysis, self.channel, "ggZZ", self.systematic, self.production),
                 TemplateFromFile(2, self.blindstatus, self.normalization, self.analysis, self.channel, "ZX", self.systematic, self.production),
                ]
    if self.blindstatus == "blind" or config.unblinddata:
        templates += [
                      TemplateFromFile(1, self.blindstatus, self.normalization, self.analysis, self.channel, "data", self.production, self.blindstatus) #blindstatus is here twice
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
                template.Scale(1/template.Integral())
            template.AddToLegend(legend)

    for axis in axes:
        hstack = ROOT.THStack(axis.name, axis.title)
        for template in templates:
            if template.color:
                hstack.Add(template.Projection(axis.index), template.hstackoption)
        hstack.Draw("nostack")
        hstack.GetXaxis().SetTitle(axis.title)
        legend.Draw()
        dir = os.path.join(config.plotsbasedir, "templateprojections", "blind" if self.blindstatus == "blind" else "fullrange", str(self.normalization), "{}_{}/{}".format(self.analysis, self.production, self.channel))
        try:
            os.makedirs(dir)
        except OSError:
            pass
        for ext in exts:
            c1.SaveAs(os.path.join(dir, "{}.{}".format(axis.name, ext)))

def projections(*args):
    Projections(*args).projections()

if __name__ == "__main__":
  for production in productions:
    print production
    for channel in channels:
      for analysis in analyses:
        for normalization in normalizations:
          for blindstatus in blindstatuses:
#            if channel != "2e2mu" or analysis != "fa3" or normalization != "areanormalize": continue
            projections(channel, analysis, normalization, production, blindstatus)
