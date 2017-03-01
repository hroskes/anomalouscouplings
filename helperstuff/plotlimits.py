import array
from collections import namedtuple
from combinehelpers import Luminosity
import config
from enums import Analysis, EnumItem, MyEnum
from extendedcounter import ExtendedCounter
from itertools import islice
import ROOT
import stylefunctions as style
import sys

filenametemplate = "higgsCombine_{append}{scanrangeappend}.MultiDimFit.mH125.root"
branchname = "CMS_zz4l_fai1"

Scan = namedtuple("Scan", "name title color style")

def plotlimits(outputfilename, analysis, *args, **kwargs):
    analysis = Analysis(analysis)
    productions = None
    legendposition = (.2, .7, .6, .9)
    moreappend = ""
    luminosity = None
    scanranges = None
    nuisance = None
    yaxistitle = "-2#Deltaln L"
    for kw, kwarg in kwargs.iteritems():
        if kw == "productions":
            productions = kwarg
        elif kw == "legendposition":
            legendposition = kwarg
        elif kw == "CLtextposition":
            CLtextposition = kwarg
        elif kw == "moreappend":
            moreappend = kwarg
        elif kw == "luminosity":
            luminosity = kwarg
        elif kw == "scanranges":
            scanranges = kwarg
        elif kw == "nuisance":
            nuisance = kwarg
        elif kw == "yaxistitle":
            yaxistitle = kwarg
        else:
            raise TypeError("Unknown kwarg {}={}".format(kw, kwarg))

    xmin = min(scanrange[1] for scanrange in scanranges)
    xmax = max(scanrange[2] for scanrange in scanranges)

    scans = []
    uptocolor = 1
    for arg in args:
        if arg == "obs":
            scans.append(Scan("obs{}".format(moreappend), "Observed, {} = 0 or #pi".format(analysis.phi), 1, 1))
            if productions is None:
                raise ValueError("No productions provided!")
        else:
            try:
                arg = float(arg)
            except ValueError:
                raise TypeError("Extra arguments to plotlimits have to be 'obs' or a float!")
            if arg == 0:
                scans.append(Scan("exp_{}{}".format(arg, moreappend), "Expected, {} = 0 or #pi".format(analysis.phi), uptocolor, 2))
            else:
                scans.append(Scan("exp_{}{}".format(arg, moreappend), "Expected, {} = {:+.2f}, {} = 0 or #pi".format(analysis.title(), arg, analysis.phi).replace("+", "#plus ").replace("-", "#minus "), uptocolor, 2))
            uptocolor += 1

    if luminosity is None:
        if productions is None or not config.unblindscans:
            luminosity = float(Luminosity("forexpectedscan"))
        else:
            luminosity = sum(float(Luminosity("fordata", production)) for production in productions)

    mg = ROOT.TMultiGraph()
    l = ROOT.TLegend(*legendposition)
    l.SetFillStyle(0)
    l.SetBorderSize(0)

    for scan in scans:
        NLL = ExtendedCounter()
        for scanrange in scanranges:
            if scanrange == (100, -1, 1):
                scanrangeappend = ""
            else:
                scanrangeappend = "_{},{},{}".format(*scanrange)
            f = ROOT.TFile(filenametemplate.format(append=scan.name, scanrangeappend=scanrangeappend))
            assert f
            t = f.Get("limit")
            assert t


            for entry in islice(t, 1, None):
                fa3 = getattr(t, branchname)
                if nuisance is None:
                    deltaNLL = t.deltaNLL+t.nll
                    NLL[fa3] = 2*deltaNLL
                else:
                     NLL[fa3] = getattr(t, nuisance)
            if 1 not in NLL and -1 in NLL:
                NLL[1] = NLL[-1]

        c1 = ROOT.TCanvas("c1", "", 8, 30, 800, 800)
        if nuisance is None: NLL.zero()
        g = NLL.TGraph()
        mg.Add(g)

        g.SetLineColor(scan.color)
        g.SetLineStyle(scan.style)
        g.SetLineWidth(2)

        l.AddEntry(g, scan.title, "l")

    mg.Draw("AL")
    mg.GetXaxis().SetTitle("{} cos({})".format(analysis.title(), analysis.phi_lower))
    mg.GetXaxis().SetRangeUser(-1, 1)
    mg.GetYaxis().SetTitle(yaxistitle)
    l.Draw()

    style.applycanvasstyle(c1)
    style.applyaxesstyle(mg)
    style.CMS("Preliminary", luminosity)

    if nuisance is None:
        drawlines(CLtextposition, xmin=xmin, xmax=xmax)
    for ext in "png eps root pdf".split():
        outputfilename = outputfilename.replace("."+ext, "")
    for ext in "png eps root pdf".split():
        c1.SaveAs("{}.{}".format(outputfilename, ext))

#https://root.cern.ch/phpBB3/viewtopic.php?t=10159
def GetNDC(x, y):
    ROOT.gPad.Update()#this is necessary!
    xndc = (x - ROOT.gPad.GetX1())/(ROOT.gPad.GetX2()-ROOT.gPad.GetX1())
    yndc = (y - ROOT.gPad.GetY1())/(ROOT.gPad.GetY2()-ROOT.gPad.GetY1())
    return xndc, yndc

cache = []

class XPos(MyEnum):
    enumitems = (
                 EnumItem("right"),
                 EnumItem("left"),
                 EnumItem("custom"),
                )
    def __init__(self, value):
        try:
            value = float(value)
        except (TypeError, ValueError):
            super(XPos, self).__init__(value)
        else:
            super(XPos, self).__init__("custom")
            self.custompos = value
    def __nonzero__(self):
        return self == "right" or self == "left" or -1 <= self.custompos <= 1
    def TPaveText(self, ypos):
        xsize = .1
        ysize = .03
        yshift = 0
        if self == "right":
            x2, y1 = GetNDC(1, ypos)
            x1 = x2 - xsize
        elif self == "left":
            x1, y1 = GetNDC(-1, ypos)
            x2 = x1 + xsize
        elif self == "custom":
            x1, y1 = GetNDC(self.custompos, ypos)
            x2 = x1 + xsize
        else:
            assert False

        y1 += yshift  #make some room between the text and the line
        y2 = y1 + ysize

        ymax = 1-ROOT.gPad.GetTopMargin()
        if y1 < ymax < y2:
            y1 -= (y2-ymax)
            y2 -= (y2-ymax)
        if y1 > ymax:
            y1 += 100  #so it doesn't get in the CMS header
            y2 += 100

        return ROOT.TPaveText(x1, y1, x2, y2, "NDC")

def drawlines(xpostext="left", xmin=-1, xmax=1):
    xpostext = XPos(xpostext)

    line68 = ROOT.TLine()
    line68.SetLineStyle(9)
    line68.DrawLine(xmin,1,xmax,1)
    cache.append(line68)
    line95 = ROOT.TLine()
    line95.SetLineStyle(9)
    line95.DrawLine(xmin,3.84,xmax,3.84)
    cache.append(line95)

    oneSig = xpostext.TPaveText(1)
    oneSig.SetFillColor(0)
    oneSig.SetFillStyle(0)
    oneSig.SetTextFont(42)
    oneSig.SetBorderSize(0)

    twoSig = xpostext.TPaveText(3.84)
    twoSig.SetFillColor(0)
    twoSig.SetFillStyle(0)
    twoSig.SetTextFont(42)
    twoSig.SetBorderSize(0)

    if xpostext:
        oneSig.AddText("68% CL")
        twoSig.AddText("95% CL")

    oneSig.Draw()
    cache.append(oneSig)
    twoSig.Draw()
    cache.append(twoSig)
