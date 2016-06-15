import collections
import helperstuff.style
from helperstuff.enums import channels
from printrates import tfiles
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
    def __init__(self, name, title, color, infile, isbkg):
        self.name = name
        self.title = title
        self.color = color
        self.infile = infile
        self.isbkg = isbkg
        self.h = None
        self.projections = {}

    def GetFromFile(self, channel):
        if self.infile:
            self.h = getattr(tfiles[channel.templatesfile(self.isbkg)], self.name)

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

exts = "png", "eps", "root", "pdf"

def projections(channel):
    templates = [
                 Template("template0PlusAdapSmoothMirror", "0^{+}", 1, True, False),
                 Template("template0MinusAdapSmoothMirror", "0^{-}", ROOT.kCyan, True, False),
                 Template("templateIntAdapSmoothMirror", "", 0, True, False),
                 Template("templateMix", "f_{a3}=0.5", 3, False, False),
                 Template("templateMixMinus", "f_{a3}=-0.5", 4, False, False),
                 Template("templateqqZZAdapSmoothMirror", "qq#rightarrowZZ", 6, True, True),
                 Template("templateggZZAdapSmoothMirror", "gg#rightarrowZZ", ROOT.kOrange+6, True, True),
                 Template("templateZXAdapSmoothMirror", "Z+X", 2, True, True),
                ]

    for i, template in enumerate(templates):
        template.GetFromFile(channel)

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
            hstack.Add(template.Projection(axis.index))
            template.AddToLegend(legend)
        hstack.Draw("nostack hist")
        hstack.GetXaxis().SetTitle(axis.title)
        legend.Draw()
        for ext in exts:
            c1.SaveAs("~/www/anomalouscouplings/templateprojections/{}/{}.{}".format(channel, axis.name, ext))

if __name__ == "__main__":
  for channel in channels:
    projections(channel)
