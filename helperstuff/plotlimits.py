#!/usr/bin/env python

import array
import math
import os
import sys

from collections import namedtuple
from itertools import islice

import ROOT

import config
import stylefunctions as style

from combinehelpers import Luminosity
from enums import Analysis, EnumItem, MyEnum
from extendedcounter import ExtendedCounter
from utilities import cache

filenametemplate = "higgsCombine_{append}{scanrangeappend}.MultiDimFit.mH125.root"

Scan = namedtuple("Scan", "name title color style")

def plottitle(nuisance):
    nuisance = actualvariable(nuisance)
    if nuisance == "CMS_zz4l_fai1": return "fai"
    if nuisance == "muV_scaled": return "muV"
    if nuisance == "muf_scaled": return "muf"
    if nuisance == "JES": return nuisance
    assert False, nuisance
def xaxistitle(POI, analysis):
    POI = actualvariable(POI)
    if POI == "CMS_zz4l_fai1": return "{} cos({})".format(analysis.title(), analysis.phi_lower)
    if POI == "muV_scaled": return "#mu_{V}"
    if POI == "muf_scaled": return "#mu_{f}"
    if nuisance == "JES": return nuisance
    assert False
def xaxisrange(POI):
    POI = actualvariable(POI)
    if POI == "CMS_zz4l_fai1": return -1.0, 1.0
def yaxistitle(nuisance, analysis):
    nuisance = actualvariable(nuisance)
    if nuisance is None: return "#minus2 #Deltaln L"
    return xaxistitle(POI=nuisance, analysis=analysis)
def actualvariable(variable):
    if variable in ("r_VVH", "muV", "RV"): return "muV_scaled"
    if variable in ("r_ffH", "muf", "RF"): return "muf_scaled"
    return variable

def plotlimits(outputfilename, analysis, *args, **kwargs):
    if not args: return
    analysis = Analysis(analysis)
    productions = None
    legendposition = (.2, .7, .6, .9)
    moreappend = ""
    luminosity = None
    scanranges = None
    nuisance = None
    POI = "CMS_zz4l_fai1"
    fixfai = False
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
            nuisance = actualvariable(kwarg)
        elif kw == "POI":
            POI = actualvariable(kwarg)
        elif kw == "fixfai":
            fixfai = kwarg
        else:
            raise TypeError("Unknown kwarg {}={}".format(kw, kwarg))

    xmin = min(scanrange[1] for scanrange in scanranges)
    xmax = max(scanrange[2] for scanrange in scanranges)

    scans = []
    uptocolor = 1
    for arg in args:
        print fixfai
        if fixfai:
            phipart = "{} = 0".format(analysis.title())
        else:
            phipart = "{} = 0 or #pi".format(analysis.phi)
        if arg == "obs":
            scans.append(Scan("obs{}".format(moreappend), "Observed, {}".format(phipart), 1, 1))
            print "Observed, {}".format(phipart)
            if productions is None:
                raise ValueError("No productions provided!")
        else:
            try:
                arg = float(arg)
            except ValueError:
                raise TypeError("Extra arguments to plotlimits have to be 'obs' or a float!")
            if arg == 0:
                scans.append(Scan("exp_{}{}".format(arg, moreappend), "Expected, {}".format(phipart, uptocolor, 2)))
            else:
                assert not fixfai
                scans.append(Scan("exp_{}{}".format(arg, moreappend), "Expected, {} = {:+.2f}, {}".format(analysis.title(), arg, phipart).replace("+", "#plus ").replace("-", "#minus "), uptocolor, 2))
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
                fa3 = getattr(t, POI)
                if nuisance is None:
                    deltaNLL = t.deltaNLL+t.nll+t.nll0
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
    mg.GetXaxis().SetTitle(xaxistitle(POI, analysis))
    if xaxisrange(POI) is not None:
        mg.GetXaxis().SetRangeUser(*xaxisrange(POI))
    mg.GetYaxis().SetTitle(yaxistitle(nuisance, analysis))
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
def GetNDC(x, y, logscale=False):
    ROOT.gPad.Update()#this is necessary!
    if logscale:
        y = math.log10(y)

    xndc = (x - ROOT.gPad.GetX1())/(ROOT.gPad.GetX2()-ROOT.gPad.GetX1())
    yndc = (y - ROOT.gPad.GetY1()) / (ROOT.gPad.GetY2() - ROOT.gPad.GetY1())
    return xndc, yndc

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
    def TPaveText(self, ypos, logscale=False, yshift=.01, xsize=.1, ysize=.03):
        xsize = .1
        ysize = .03
        if self == "right":
            x2, y1 = GetNDC(.95, ypos, logscale)
            x1 = x2 - xsize
        elif self == "left":
            x1, y1 = GetNDC(-.95, ypos, logscale)
            x2 = x1 + xsize
        elif self == "custom":
            x1, y1 = GetNDC(self.custompos, ypos, logscale)
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

@cache
def drawlines(xpostext="left", xmin=-1, xmax=1, logscale=False, xsize=.1, ysize=.03, yshift=0, yshift68=None, yshift95=None, textsize=None, PRL=False, CLbelow=False, arbitraryparameter=None):
    """
    xpostext: "left", "right", or a float, determines where the text 68% CL and 95% CL goes
    xmin: minimum of plot
    xmax: maximum of plot
    logscale: is the plot in logscale
    xsize, ysize, textsize: xsize of text
    yshift, yshift68, yshift95: y shift of text from lines
    PRL: quick switch for big text for PRL
      CLbelow: for PRL, put CL text below the lines rather than above
    arbitraryparameter: doesn't do anything, but the lines are stored separately in the cache
    """
    xpostext = XPos(xpostext)
    line68 = ROOT.TLine()
    line68.SetLineStyle(9)
    line68.DrawLine(xmin,1,xmax,1)
    line95 = ROOT.TLine()
    line95.SetLineStyle(9)
    line95.DrawLine(xmin,3.84,xmax,3.84)

    if CLbelow and not PRL:
        raise ValueError("CLbelow only works with PRL.  Turn it on, or set yshift yourself")
    if PRL and not (xsize == .1 and ysize == .03 and yshift == 0 and textsize is yshift68 is yshift95 is None):
        raise ValueError("To set PRL options, keep xsize, ysize, yshift, and textsize as their defaults!")

    if PRL:
        xsize = .15
        textsize = .044
        if CLbelow:
            yshift=-.04
        else:
            yshift=.01

    if yshift68 is None: yshift68 = yshift
    if yshift95 is None: yshift95 = yshift

    oneSig = xpostext.TPaveText(1, logscale=logscale, xsize=xsize, ysize=ysize, yshift=yshift68)
    oneSig.SetFillColor(0)
    oneSig.SetFillStyle(0)
    oneSig.SetTextFont(42)
    oneSig.SetBorderSize(0)

    twoSig = xpostext.TPaveText(3.84, logscale=logscale, xsize=xsize, ysize=ysize, yshift=yshift95)
    twoSig.SetFillColor(0)
    twoSig.SetFillStyle(0)
    twoSig.SetTextFont(42)
    twoSig.SetBorderSize(0)

    if textsize is not None:
        oneSig.SetTextSize(textsize)
        twoSig.SetTextSize(textsize)

    if xpostext:
        oneSig.AddText("68% CL")
        twoSig.AddText("95% CL")

    oneSig.Draw()
    twoSig.Draw()

    return line68, line95, oneSig, twoSig

if __name__ == "__main__":
    args, kwargs = [], {}
    for arg in sys.argv[1:]:
        if "=" in arg:
            kw, kwarg = arg.split("=")
            if kw in kwargs:
                raise TypeError("Duplicate kwarg {}!".format(kw))
            if kw == "scanranges":
                kwarg = kwarg.replace(";", ":")
                kwarg = [tuple(float(_2) for _2 in _.split(",")) for _ in kwarg.split(":")]
                assert all(len(_) == 3 for _ in kwarg)
            if kw == "productions":
                kwarg = kwarg.split(",")
            kwargs[kw] = kwarg
        else:
            args.append(arg)

    plotlimits(*args, **kwargs)

    outputfilename = kwargs.get("outputfilename", args[0])

    for ext in "png eps root pdf".split():
        outputfilename = outputfilename.replace("."+ext, "")
    with open(outputfilename+".txt", 'w') as f:
        f.write("cd {} &&\n".format(os.getcwd()))
        f.write("python " + " ".join(sys.argv))

@cache
def arrowatminimum(graph, abovexaxis=True):
  minimum = (float("nan"), float("inf"))
  for i, x, y in zip(xrange(graph.GetN()), graph.GetX(), graph.GetY()):
    if y < minimum[1]:
      minimum = x, y

  x = minimum[0]
  if abovexaxis:
    y1, y2 = .15, .1
  else:
    y1, y2 = .1**2/.15, .1
  arrow = ROOT.TArrow(x, y1, x, y2, .01, "|>")
  arrow.SetLineStyle(graph.GetLineStyle())
  arrow.SetLineColor(graph.GetLineColor())
  arrow.SetFillColor(graph.GetLineColor())
  arrow.SetLineWidth(graph.GetLineWidth())
  return arrow
