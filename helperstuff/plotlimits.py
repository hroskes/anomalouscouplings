#!/usr/bin/env python

import array
import math
import os
import random
import re
import sys

from collections import namedtuple
from itertools import islice, izip

import numpy as np

import ROOT

import config
import stylefunctions as style

from combinehelpers import Luminosity, mixturesign, sigmaioversigma1
from enums import Analysis, Category, channels, EnumItem, MyEnum, ProductionMode
from samples import ReweightingSample, Sample, samplewithfai
from extendedcounter import ExtendedCounter
from samples import samplewithfeLfeR
from utilities import cache, cd
from yields import YieldValue

filenametemplate = "higgsCombine_{append}{scanrangeappend}.MultiDimFit.mH125.root"

class Scan(namedtuple("Scan", "name title color style badpoints")):
    def __new__(cls, name, title, color, style, badpoints=None):
        return super(Scan, cls).__new__(cls, name, title, color, style, badpoints=badpoints)

def plottitle(nuisance):
    nuisance = actualvariable(nuisance)
    if nuisance == "CMS_zz4l_fai1": return "fai"
    if nuisance == "CMS_zz4l_fai2": return "faj"
    if nuisance == "RV": return "muV"
    if nuisance == "RF": return "muf"
    if nuisance == "CMS_scale_j_13TeV_2016": return nuisance
    assert False, nuisance
def xaxistitle(POI, analysis, faifor=None):
    POI = actualvariable(POI)
    if faifor == "decay": faifor = None
    if POI == "CMS_zz4l_fai1": return "{} cos({})".format(analysis.title(superscript=faifor), analysis.phi_lower)
    if POI == "CMS_zz4l_fai2": return "{} cos({})".format(analysis.fais[1].title(superscript=faifor), analysis.fais[1].phi_lower)
    if POI == "CMS_zz4l_fai3": return "{} cos({})".format(analysis.fais[2].title(superscript=faifor), analysis.fais[2].phi_lower)
    if POI == "CMS_zz4l_fai4": return "{} cos({})".format(analysis.fais[3].title(superscript=faifor), analysis.fais[3].phi_lower)
    if POI == "RV": return "#mu_{V}"
    if POI == "RF": return "#mu_{f}"
    if nuisance == "CMS_scale_j_13TeV_2016": return nuisance
    assert False
def xaxisrange(POI):
    POI = actualvariable(POI)
    if POI == "CMS_zz4l_fai1": return -1.0, 1.0
    if POI == "CMS_zz4l_fai2": return -1.0, 1.0
    if POI == "CMS_zz4l_fai3": return -1.0, 1.0
    if POI == "CMS_zz4l_fai4": return -1.0, 1.0
def yaxistitle(nuisance, analysis):
    nuisance = actualvariable(nuisance)
    if nuisance is None: return "#minus2 #Deltaln L"
    return xaxistitle(POI=nuisance, analysis=analysis)
def actualvariable(variable):
    #if variable in ("r_VVH", "muV", "RV"): return "muV_scaled"
    #if variable in ("r_ffH", "muf", "RF"): return "muf_scaled"
    return variable


@cache
def sigmai_VBFreco(analysis, production):
    #copied from step12_categorysystematics/latextable.py
    t = ROOT.TChain("candTree")
    for sample in ProductionMode("VBF").allsamples(production):
        t.Add(sample.withdiscriminantsfile())
    weightparts = [
        "ZZMass>{}".format(config.m4lmin),
        "ZZMass<{}".format(config.m4lmax),
        ReweightingSample("VBF", analysis.purehypotheses[1]).weightname(),
        "||".join("(category_{}=={})".format(analysis.categoryname, idnumber) for idnumber in Category("VBFtagged").idnumbers),
    ]
    wt = " * ".join("("+_+")" for _ in weightparts)
    hname = "h{}".format(random.getrandbits(100))
    success = t.Draw("1>>{}".format(hname), wt)
    assert success
    return getattr(ROOT, hname).Integral() / sum(
        Sample.effectiveentries(
            reweightfrom=reweightfrom,
            reweightto=ReweightingSample("VBF", analysis.purehypotheses[1]),
        ) for reweightfrom in ProductionMode("VBF").allsamples(production)
    ) * ReweightingSample("VBF", analysis.purehypotheses[1]).xsec / ReweightingSample("VBF", analysis.purehypotheses[0]).xsec

@cache
def fai_VBFreco(analysis, production, faidecay):
    s = samplewithfai("ggH", analysis, faidecay)
    a1 = s.g1
    ai = getattr(s, analysis.couplingname)
    sigma1 = sum(YieldValue("VBF", analysis, "VBFtagged", c).value for c in channels)
    sigmai = sigmai_VBFreco(analysis, production)
    sigmai /= sigmaioversigma1(analysis, "ggH")
    print sigma1, sigmai, a1, ai
    return mixturesign(analysis) * math.copysign(1, a1*ai) * ai**2*sigmai / (a1**2*sigma1 + ai**2*sigmai)



def plotlimits(outputfilename, analysis, *args, **kwargs):
    if not args: return
    analysis = Analysis(analysis)
    productions = None
    legendposition = (.2, .7, .6, .9)
    CLtextposition = "left"
    moreappend = ""
    luminosity = None
    scanranges = None
    nuisance = None
    POI = "CMS_zz4l_fai1"
    fixfai = False
    CMStext = "Preliminary"
    drawCMS = True
    contactterms = False
    infilename = filenametemplate
    xtitle = None
    scanfai = analysis
    killpoints = []
    faifor = "decay"
    xaxislimits = None
    adddirectories = []
    scale = 1
    useNLLandNLL0 = True
    plotcopier = ROOT
    combinelogs = None
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
        elif kw == "CMStext":
            CMStext = kwarg
        elif kw == "drawCMS":
            drawCMS = kwarg
        elif kw == "scanfai":
            scanfai = kwarg
        elif kw == "contactterms":
            contactterms = kwarg
        elif kw == "infilename":
            infilename = kwarg
        elif kw == "xtitle":
            xtitle = kwarg
        elif kw == "killpoints":
            if isinstance(kwarg, basestring):
                kwarg = kwarg.replace(":", ";")
                kwarg = kwarg.split(";")
            for _ in kwarg:
                if isinstance(_, basestring): _ = [float(a) for a in _.split(",")]
                assert len(_) == 2
                killpoints.append(_)
        elif kw == "faifor":
            faifor = kwarg
        elif kw == "xaxislimits":
            xaxislimits = kwarg
        elif kw == "adddirectories":
            adddirectories = kwarg
        elif kw == "scale":
            scale = kwarg
        elif kw == "useNLLandNLL0":
            useNLLandNLL0 = bool(int(kwarg))
        elif kw == "plotcopier":
            plotcopier = kwarg
        elif kw == "combinelogs":
            combinelogs = kwarg
        else:
            raise TypeError("Unknown kwarg {}={}".format(kw, kwarg))

    if scanfai == analysis.fais[0] or analysis.dimensions == 1: pass
    elif scanfai == analysis.fais[1]: POI = "CMS_zz4l_fai2"
    elif scanfai == analysis.fais[2]: POI = "CMS_zz4l_fai3"
    elif scanfai == analysis.fais[3]: POI = "CMS_zz4l_fai4"
    else: assert False

    xmin = min(scanrange[1] for scanrange in scanranges)
    xmax = max(scanrange[2] for scanrange in scanranges)

    scans = []
    uptocolor = 1

    badpoints = [[] * len(args)]
    if combinelogs is not None:
        assert len(combinelogs) == len(args)
        badpoints = [findbadpoints(*combineloglist) for combineloglist in combinelogs]

    for arg, badpts in izip(args, badpoints):
        if fixfai:
            phipart = "{} = 0".format(analysis.title())
        else:
            phipart = "{} = 0 or #pi".format(analysis.phi)
        if arg == "obs":
            scans.append(Scan("obs{}".format(moreappend), "Observed, {}".format(phipart), 1, 1, badpoints=badpts))
            if productions is None:
                raise ValueError("No productions provided!")
        else:
            try:
                arg = float(arg)
            except ValueError:
                raise TypeError("Extra arguments to plotlimits have to be 'obs' or a float!")
            if arg == 0:
                scans.append(Scan("exp_{}{}".format(arg, moreappend), "Expected, {}".format(phipart, uptocolor, 2), uptocolor, 2, badpoints=badpts))
            else:
                assert not fixfai
                scans.append(Scan("exp_{}{}".format(arg, moreappend), "Expected, {} = {:+.2f}, {}".format(analysis.title(), arg, phipart).replace("+", "#plus ").replace("-", "#minus "), uptocolor, 2, badpoints=badpts))
            uptocolor += 1

    if luminosity is None:
        luminosity = sum(float(Luminosity("fordata", production)) for production in productions)

    mg = ROOT.TMultiGraph()
    l = ROOT.TLegend(*legendposition)
    l.SetFillStyle(0)
    l.SetBorderSize(0)

    for scan in scans:
        finalNLL = ExtendedCounter()
        for directory in ["."] + adddirectories:
            usemoreappend = ""
            if directory != ".":
                scale = 1
                if isinstance(directory, tuple): directory, usemoreappend, scale = directory
            with cd(directory):
                NLL = ExtendedCounter()
                for scanrange in scanranges:
                    if scanrange == (101, -1, 1):
                        scanrangeappend = ""
                    else:
                        scanrangeappend = "_{},{},{}".format(*scanrange)
                    f = ROOT.TFile(infilename.format(append=scan.name+usemoreappend, scanrangeappend=scanrangeappend))
                    assert f
                    t = f.Get("limit")
                    if not t: print os.getcwd(); f.ls()
                    assert t


                    t.GetEntry(0)
                    if getattr(t, POI) == -1: startfrom = 0
                    else: startfrom = 1

                    for entry in islice(t, startfrom, None):
                        fa3 = getattr(t, POI)
                        if POI == "CMS_zz4l_fai1":
                            if faifor == "decay": pass
                            elif faifor in ("VBF", "VH", "ZH", "WH"): fa3 = samplewithfai("ggH", analysis, fa3).fai(faifor, analysis)
                            elif faifor in ("VBF+dec", "VH+dec", "ZH+dec", "WH+dec"): fa3 = samplewithfai("ggH", analysis, fa3).fai(faifor, analysis, withdecay=True)
                            elif faifor == "VBFreco": assert len(productions) == 1; fa3 = fai_VBFreco(analysis, productions[0], fa3)
                            else: assert False
                            if xaxislimits is not None and not (xaxislimits[0] <= fa3 <= xaxislimits[1]): continue
                        if any(_[0] < fa3 < _[1] for _ in killpoints): continue
                        if any(np.isclose(fa3, _) for _ in scan.badpoints): continue
                        if nuisance is None:
                            deltaNLL = t.deltaNLL
                            if useNLLandNLL0: deltaNLL += t.nll+t.nll0
                            if math.isinf(deltaNLL) or math.isnan(deltaNLL):
                                t.Show()
                                raise ValueError("one of the NLL variables is infinity or NaN, see ^^^^")
                            NLL[fa3] = 2*deltaNLL
                        else:
                            NLL[fa3] = getattr(t, nuisance)
                    if 1 not in NLL and -1 in NLL and xaxislimits is None and POI == "CMS_zz4l_fai1":
                        NLL[1] = NLL[-1]

                NLL *= scale

                if set(finalNLL):
                    for point in finalNLL.copy():
                        if point not in NLL:
                            try:
                                left = max(_ for _ in NLL if _<point)
                            except ValueError:
                                left = None
                            try:
                                right = min(_ for _ in NLL if _>point)
                            except ValueError:
                                right = None

                            if right is None:
                                right = left
                                left = max(_ for _ in NLL if _<right)
                            if left is None:
                                left = right
                                right = min(_ for _ in NLL if _>left)

                            NLL[point] = -(point * (NLL[right]-NLL[left]) + right * NLL[left] - left*NLL[right]) / (left - right)
                    for point in NLL.copy():
                        if point not in finalNLL:
                            del NLL[point]

            finalNLL += NLL

        c1 = plotcopier.TCanvas("c1", "", 8, 30, 800, 800)
        if nuisance is None: finalNLL.zero()
        g = finalNLL.TGraph()
        mg.Add(g)

        g.SetLineColor(scan.color)
        g.SetLineStyle(scan.style)
        g.SetLineWidth(2)

        l.AddEntry(g, scan.title, "l")

    mg.Draw("AL")
    if xtitle is None: xtitle = xaxistitle(POI, analysis, faifor=faifor)
    mg.GetXaxis().SetTitle(xtitle)
    if xaxislimits is not None:
        mg.GetXaxis().SetRangeUser(*xaxislimits)
    if xaxisrange(POI) is not None:
        mg.GetXaxis().SetRangeUser(*xaxisrange(POI))
    mg.GetYaxis().SetTitle(yaxistitle(nuisance, analysis))
    l.Draw()

    style.applycanvasstyle(c1)
    style.applyaxesstyle(mg)
    style.CMS(CMStext, luminosity, drawCMS=drawCMS)

    if nuisance is None:
        drawlines(CLtextposition, xmin=xmin, xmax=xmax)
    for ext in "png eps root pdf C".split():
        outputfilename = outputfilename.replace("."+ext, "")
    for ext in "png eps root pdf C".split():
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

def plotlimits2D(outputfilename, analysis, *args, **kwargs):
    if not args: return
    analysis = Analysis(analysis)
    productions = None
    legendposition = (.2, .7, .6, .9)
    moreappend = ""
    luminosity = None
    scanranges = None
    nuisance = None
    POI = "CMS_zz4l_fai1"
    POI2 = "CMS_zz4l_fai2"
    fixfai = False
    CMStext = "Preliminary"
    drawCMS = True
    plotcopier = ROOT
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
        elif kw == "POI2":
            POI2 = actualvariable(kwarg)
        elif kw == "fixfai":
            fixfai = kwarg
        elif kw == "CMStext":
            CMStext = kwarg
        elif kw == "drawCMS":
            drawCMS = kwarg
        elif kw == "scanfai":
            scanfai = kwarg
        elif kw == "plotcopier":
            plotcopier = kwarg
        else:
            raise TypeError("Unknown kwarg {}={}".format(kw, kwarg))

    xmin = ymin = min(scanrange[1] for scanrange in scanranges)
    xmax = ymax = max(scanrange[2] for scanrange in scanranges)

    scans = []
    uptocolor = 1
    assert len(args) == 1, args
    for arg in args:
        if fixfai:
            phipart = "{} = 0".format(analysis.title())
        else:
            phipart = "{} = 0 or #pi".format(analysis.phi)
        if arg == "obs":
            scans.append(Scan("obs{}".format(moreappend), "Observed, {}".format(phipart), 1, 1))
            if productions is None:
                raise ValueError("No productions provided!")
        else:
            try:
                arg = float(arg)
            except ValueError:
                raise TypeError("Extra arguments to plotlimits have to be 'obs' or a float!")
            if arg == 0:
                scans.append(Scan("exp_{}{}".format(arg, moreappend), "Expected, {}".format(phipart, uptocolor, 2), uptocolor, 2))
            else:
                assert not fixfai
                scans.append(Scan("exp_{}{}".format(arg, moreappend), "Expected, {} = {:+.2f}, {}".format(analysis.title(), arg, phipart).replace("+", "#plus ").replace("-", "#minus "), uptocolor, 2))
            uptocolor += 1

    if luminosity is None:
        luminosity = sum(float(Luminosity("fordata", production)) for production in productions)

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


            t.GetEntry(0)
            if getattr(t, POI) == -1: startfrom = 0
            else: startfrom = 1

            for entry in islice(t, startfrom, None):
                if abs(t.CMS_zz4l_fai1) + abs(t.CMS_zz4l_fai2) > 1+1e-14: continue
                fa3 = getattr(t, POI), getattr(t, POI2)
                if nuisance is None:
                    deltaNLL = t.deltaNLL+t.nll+t.nll0
                    NLL[fa3] = 2*deltaNLL
                else:
                     NLL[fa3] = getattr(t, nuisance)

        c1 = plotcopier.TCanvas("c1", "", 8, 30, 800, 800)
        if nuisance is None: NLL.zero()
        g = NLL.TGraph2D()

        g.SetLineColor(scan.color)
        g.SetLineStyle(scan.style)
        g.SetLineWidth(2)

        l.AddEntry(g, scan.title, "l")

    g.Draw("COLZ")
    g.GetHistogram().GetXaxis().SetTitle(xaxistitle(POI, analysis))
    if xaxisrange(POI) is not None:
        g.GetXaxis().SetRangeUser(*xaxisrange(POI))
    g.GetHistogram().GetYaxis().SetTitle(xaxistitle(POI2, analysis))
    if xaxisrange(POI2) is not None:
        g.GetYaxis().SetRangeUser(*xaxisrange(POI2))
    g.GetHistogram().GetZaxis().SetTitle(yaxistitle(nuisance, analysis))
    g.GetHistogram().GetZaxis().SetTitleOffset(1.4)
    g.GetHistogram().GetZaxis().SetRangeUser(0, 300)

    gL, gLpoint, gR, gRpoint = feLfeRgraphs()

    mg = ROOT.TMultiGraph()
    mg.Add(gL, "C")
    mg.Add(gLpoint, "P")
    mg.Add(gR, "C")
    mg.Add(gRpoint, "P")
    mg.Draw()

    style.applycanvasstyle(c1)
    c1.SetRightMargin(0.15)
    style.applyaxesstyle(g)
    style.CMS(CMStext, luminosity, drawCMS=drawCMS)

    for ext in "png eps root pdf".split():
        outputfilename = outputfilename.replace("."+ext, "")
    for ext in "png eps root pdf".split():
        c1.SaveAs("{}.{}".format(outputfilename, ext))

@cache
def feLfeRgraphs():
    fL1 = array.array("d", (samplewithfeLfeR(i/100., 0).fai("ggH", "fL1") for i in range(-100, 101)))
    fL1Zg = array.array("d", (samplewithfeLfeR(i/100., 0).fai("ggH", "fL1Zg") for i in range(-100, 101)))
    gL = ROOT.TGraph(len(fL1), fL1, fL1Zg)
    gL.SetLineColor(2)
    gL.SetLineWidth(4)

    fL1 = array.array("d", (samplewithfeLfeR(i/100., 0).fai("ggH", "fL1") for i in range(100, 101)))
    fL1Zg = array.array("d", (samplewithfeLfeR(i/100., 0).fai("ggH", "fL1Zg") for i in range(100, 101)))
    gLpoint = ROOT.TGraph(len(fL1), fL1, fL1Zg)
    gLpoint.SetMarkerColor(2)
    gLpoint.SetMarkerStyle(20)
    gLpoint.SetMarkerSize(gLpoint.GetMarkerSize()+2)

    fL1 = array.array("d", (samplewithfeLfeR(0, i/100.).fai("ggH", "fL1") for i in range(-100, 101)))
    fL1Zg = array.array("d", (samplewithfeLfeR(0, i/100.).fai("ggH", "fL1Zg") for i in range(-100, 101)))
    for _ in zip(range(-100, 101), fL1, fL1Zg): print _
    gR = ROOT.TGraph(len(fL1), fL1, fL1Zg)
    gR.SetLineColor(3)
    gR.SetLineWidth(4)

    fL1 = array.array("d", (samplewithfeLfeR(0, i/100.).fai("ggH", "fL1") for i in range(100, 101)))
    fL1Zg = array.array("d", (samplewithfeLfeR(0, i/100.).fai("ggH", "fL1Zg") for i in range(100, 101)))
    gRpoint = ROOT.TGraph(len(fL1), fL1, fL1Zg)
    gRpoint.SetMarkerColor(3)
    gRpoint.SetMarkerStyle(20)
    gRpoint.SetMarkerSize(gRpoint.GetMarkerSize()+2)

    return gL, gLpoint, gR, gRpoint

def findbadpoints(*combinelogs):
    result = set()
    for combinelog in combinelogs:
        currentpoint = None
        with open(combinelog) as f:
            for line in f:
                #https://stackoverflow.com/a/4703508/5228524
                match = re.search(r"Point [0-9]+/[0-9]+ [A-Za-z0-9_]+ = ([-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?)", line)
                if match:
                    currentpoint = float(match.group(1))
                if "VariableMetricBuilder: matrix not pos.def." in line:
                    result.add(currentpoint)
    return result
