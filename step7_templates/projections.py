import collections
from helperstuff import config, constants
from helperstuff.combinehelpers import getrate
from helperstuff.enums import analyses, Analysis, Channel, channels, EnumItem, MultiEnum, MyEnum, Systematic, Template, TemplatesFile
from helperstuff.filemanager import tfiles
import helperstuff.style
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
                 EnumItem(""),
                 EnumItem("rescalemixtures"),
                 EnumItem("areanormalize"),
                ]
normalizations = [Normalization(normalization) for normalization in Normalization.enumitems]
    
class TemplateFromFile(TemplateForProjection, MultiEnum):
    enums = [Template, Normalization]
    def __init__(self, color, *args):
        super(TemplateFromFile, self).__init__(*args)
        self.title = self.template.title()
        self.color = color
        self.h = self.template.gettemplate().Clone()
        self.projections = {}
        if self.template.productionmode == "ggH":
            scalefactor = getrate("2e2mu", self.template.productionmode) / Template(self.template.templatesfile.analysis, "2e2mu", self.template.productionmode, "0+").gettemplate().Integral()
        else:
            scalefactor = getrate(self.template.templatesfile.channel, self.template.productionmode) / self.Integral()
        self.Scale(scalefactor)
    def check(self, *args):
        if self.normalization is None:
            self.normalization = Normalization("")
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
        analyses = set(t[0].template.templatesfile.analysis for t in templatesandfactors)
        assert len(analyses) == 1
        analyses = list(analyses)[0]
        assert abs(mixturesign) == 1
        if normalization != "rescalemixtures":
            mixturesign = 0
        self.Scale(constants.JHUXS2L2la1 / (2*constants.JHUXS2L2la1 + mixturesign*analysis.interfxsec()))

exts = "png", "eps", "root", "pdf"

def projections(channel, analysis, normalization = "", systematic = ""):
    channel = Channel(channel)
    analysis = Analysis(analysis)
    normalization = Normalization(normalization)
    systematic = Systematic(systematic)
    templates = [
                 TemplateFromFile(1, normalization, analysis.signaltemplates(channel, systematic)[0]),
                 TemplateFromFile(ROOT.kCyan, normalization, analysis.signaltemplates(channel, systematic)[1]),
                 TemplateFromFile(0, normalization, analysis.signaltemplates(channel, systematic)[2]),
                ]
    templates+= [
                 TemplateSum("ggH {}=0.5".format(analysis.title()), 3, 1, (templates[0], 1), (templates[1], 1), (templates[2], 1)),
                 TemplateSum("ggH {}=-0.5".format(analysis.title()), 4, -1, (templates[0], 1), (templates[1], 1), (templates[2], -1)),
                 TemplateFromFile(6, normalization, analysis, channel, "qqZZ", systematic),
                 TemplateFromFile(ROOT.kOrange+6, normalization, analysis, channel, "ggZZ", systematic),
                 TemplateFromFile(2, normalization, analysis, channel, "ZX", systematic),
                ]
    axes[0] = Axis(analysis.purediscriminant(), analysis.purediscriminant(title=True), 0)
    axes[1] = Axis(analysis.mixdiscriminant(), analysis.mixdiscriminant(title=True), 1)

    c1 = ROOT.TCanvas()
    legend = ROOT.TLegend(.75, .65, .9, .9)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    for template in templates:
        if template.color:
            if normalization == "areanormalize":
                template.Scale(1/template.Integral())
            template.AddToLegend(legend)

    for axis in axes:
        hstack = ROOT.THStack(axis.name, axis.title)
        for template in templates:
            if template.color:
                hstack.Add(template.Projection(axis.index))
        hstack.Draw("nostack hist")
        hstack.GetXaxis().SetTitle(axis.title)
        legend.Draw()
        dir = os.path.join(config.plotsbasedir, "templateprojections", str(normalization))
        try:
            os.makedirs(os.path.join(dir, "{}/{}".format(analysis, channel)))
        except OSError:
            pass
        for ext in exts:
            c1.SaveAs(os.path.join(dir, "{}/{}/{}.{}".format(analysis, channel, axis.name, ext)))

if __name__ == "__main__":
  for channel in channels:
    for analysis in analyses:
      for normalization in normalizations:
#        if channel != "2e2mu" or analysis != "fa3" or normalization != "areanormalize": continue
        projections(channel, analysis, normalization)
