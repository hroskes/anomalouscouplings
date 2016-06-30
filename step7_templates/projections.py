import collections
from helperstuff import config
import helperstuff.style
from helperstuff.enums import analyses, Analysis, Channel, channels, Template, TemplatesFile
from helperstuff.filemanager import tfiles
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
    
class TemplateFromFile(TemplateForProjection):
    def __init__(self, color, *args):
        self.template = Template(*args)
        self.title = self.template.title()
        self.color = color
        self.h = self.template.gettemplate()
        self.projections = {}

class TemplateSum(TemplateForProjection):
    def __init__(self, title, color, *templatesandfactors):
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

exts = "png", "eps", "root", "pdf"

def projections(channel, analysis, areanormalize=False, systematic = ""):
    channel = Channel(channel)
    analysis = Analysis(analysis)
    templates = [
                 TemplateFromFile(1, analysis.signaltemplates(channel)[0]),
                 TemplateFromFile(ROOT.kCyan, analysis.signaltemplates(channel)[1]),
                 TemplateFromFile(0, analysis.signaltemplates(channel)[2]),
                ]
    templates+= [
                 TemplateSum("ggH {}=0.5".format(analysis.title()), 3, (templates[0], 1), (templates[1], 1), (templates[2], 1)),
                 TemplateSum("ggH {}=-0.5".format(analysis.title()), 4, (templates[0], 1), (templates[1], 1), (templates[2], -1)),
                 TemplateFromFile(6, analysis, channel, "qqZZ"),
                 TemplateFromFile(ROOT.kOrange+6, analysis, channel, "ggZZ"),
                 TemplateFromFile(2, analysis, channel, "ZX"),
                ]
    axes[0] = Axis(analysis.purediscriminant(), analysis.purediscriminant(title=True), 0)
    axes[1] = Axis(analysis.mixdiscriminant(), analysis.mixdiscriminant(title=True), 1)

    c1 = ROOT.TCanvas()
    for axis in axes:
        hstack = ROOT.THStack(axis.name, axis.title)
        legend = ROOT.TLegend(.75, .65, .9, .9)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        for template in templates:
            if not template.color: continue
            if areanormalize:
                template.Scale(1 / template.Integral())
            hstack.Add(template.Projection(axis.index))
            template.AddToLegend(legend)
        hstack.Draw("nostack hist")
        hstack.GetXaxis().SetTitle(axis.title)
        legend.Draw()
        dir = os.path.join(config.plotsbasedir, "templateprojections")
        if areanormalize:
            dir = os.path.join(dir, "areanormalized")
        try:
            os.makedirs(os.path.join(dir, str(channel)))
        except OSError:
            pass
        try:
            os.makedirs(os.path.join(dir, "{}/{}".format(analysis, channel)))
        except OSError:
            pass
        for ext in exts:
            c1.SaveAs(os.path.join(dir, "{}/{}/{}.{}".format(analysis, channel, axis.name, ext)))

if __name__ == "__main__":
  for channel in channels:
    for analysis in analyses:
      projections(channel, analysis, areanormalize=False)
      projections(channel, analysis, areanormalize=True)
