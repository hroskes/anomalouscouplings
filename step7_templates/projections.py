import collections
from helperstuff import config
import helperstuff.style
from helperstuff.enums import Channel, channels
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
        Axis("D0minus", "D_{0^{-}}", 0),
        Axis("DCP", "D_{CP}", 1),
        Axis("Dbkg", "D_{bkg}", 2),
       ]

class Template(object):
    def __init__(self, name, title, color, infile, isbkg, run1name=None):
        self.name = name
        self.run1name = run1name
        if run1name is None:
            self.run1name = name
        self.title = title
        self.color = color
        self.infile = infile
        self.isbkg = isbkg
        self.h = None
        self.projections = {}

    def GetFromFile(self, channel, run1=False):
        if self.infile:
            self.h = getattr(tfiles[channel.templatesfile(self.isbkg, run1=run1)], self.name if not run1 else self.run1name)

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

exts = "png", "eps", "root", "pdf"

def projections(channel, run1=False, areanormalize=False):
    channel = Channel(channel)
    templates = [
                 Template("template0PlusAdapSmoothMirror", "0^{+}", 1, True, False),
                 Template("template0MinusAdapSmoothMirror", "0^{-}", ROOT.kCyan, True, False),
                 Template("templateIntAdapSmoothMirror", "", 0, True, False),
                 Template("templateMix", "f_{a3}=0.5", 3, False, False),
                 Template("templateMixMinus", "f_{a3}=-0.5", 4, False, False),
                 Template("templateqqZZAdapSmoothMirror", "qq#rightarrowZZ", 6, True, True, run1name="template_qqZZ"),
                 Template("templateggZZAdapSmoothMirror", "gg#rightarrowZZ", ROOT.kOrange+6, True, True, run1name="template_ggZZ"),
                 Template("templateZXAdapSmoothMirror", "Z+X", 2, True, True, run1name="template_ZX"),
                ]

    for i, template in enumerate(templates):
        template.GetFromFile(channel, run1=run1)

    if run1:
        integralSM = templates[0].Integral()
        yields = {
                  Channel("2e2mu"): 7.9085,
                  Channel("4e"): 3.1441,
                  Channel("4mu"): 6.0802,
                 }
        for template in templates[0:3]:
            template.Scale(yields[channel] / integralSM)

        #qqZZ
        yields = {
                  Channel("2e2mu"): 8.8585,
                  Channel("4e"): 2.9364,
                  Channel("4mu"): 7.6478,
                 }
        templates[5].Scale(yields[channel])

        #ggZZ
        yields = {
                  Channel("2e2mu"): 0.5005,
                  Channel("4e"): 0.2041,
                  Channel("4mu"): 0.4131,
                 }
        templates[6].Scale(yields[channel])

        #Z+X
        yields = {
                  Channel("2e2mu"): 4.2929,
                  Channel("4e"): 2.7676,
                  Channel("4mu"): 1.1878,
                 }
        templates[7].Scale(yields[channel])

    templates[3].h = templates[0].h.Clone()
    templates[3].h.Add(templates[1].h)
    templates[3].h.Add(templates[2].h)
    templates[3].Scale(.5)

    templates[4].h = templates[0].h.Clone()
    templates[4].h.Add(templates[1].h)
    templates[4].h.Add(templates[2].h, -1)
    templates[4].Scale(.5)

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
        if run1:
            dir = os.path.join(dir, "run1")
        try:
            os.makedirs(os.path.join(dir, str(channel)))
        except OSError:
            pass
        for ext in exts:
            c1.SaveAs(os.path.join(dir, "{}/{}.{}".format(channel, axis.name, ext)))

if __name__ == "__main__":
  for channel in channels:
    projections(channel, run1=False, areanormalize=False)
    projections(channel, run1=True, areanormalize=False)
    projections(channel, run1=False, areanormalize=True)
    projections(channel, run1=True, areanormalize=True)
