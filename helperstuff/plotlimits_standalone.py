#!/usr/bin/env python

import array
import math
import os
import sys

from collections import namedtuple
from itertools import islice

import ROOT

from Alignment.OfflineValidation.TkAlAllInOneTool.helperFunctions import parsecolor, parsestyle

import stylefunctions as style

from extendedcounter import ExtendedCounter
from utilities import cache

Scan = namedtuple("Scan", "filenames title color style")

def plotlimits(outputfile, *scans, **kwargs):
    legendposition = (.2, .7, .6, .9)
    CLtextposition = "left"
    luminosity = None
    scanrange = -1, 1
    POI = "CMS_zz4l_fai1"
    CMStext = "Preliminary"
    drawCMS = True
    xtitle = None
    ytitle = None
    for kw, kwarg in kwargs.iteritems():
        if kw == "legendposition":
            legendposition = kwarg
        elif kw == "CLtextposition":
            CLtextposition = kwarg
        elif kw == "luminosity":
            luminosity = kwarg
        elif kw == "scanrange":
            scanrange = kwarg
        elif kw == "POI":
            POI = actualvariable(kwarg)
        elif kw == "CMStext":
            CMStext = kwarg
        elif kw == "drawCMS":
            drawCMS = kwarg
        elif kw == "xtitle":
            xtitle = kwarg
        elif kw == "ytitle":
            ytitle = kwarg
        else:
            raise TypeError("Unknown kwarg {}={}".format(kw, kwarg))

    xmin, xmax = scanrange

    mg = ROOT.TMultiGraph()
    l = ROOT.TLegend(*legendposition)
    l.SetFillStyle(0)
    l.SetBorderSize(0)

    for scan in scans:
        NLL = ExtendedCounter()
        for filename in scan.filenames:
            f = ROOT.TFile(filename)
            assert f
            t = f.Get("limit")
            assert t


            t.GetEntry(0)
            if getattr(t, POI) == -1: startfrom = 0
            else: startfrom = 1

            for entry in islice(t, startfrom, None):
                fa3 = getattr(t, POI)
                deltaNLL = t.deltaNLL+t.nll+t.nll0
                NLL[fa3] = 2*deltaNLL
            if 1 not in NLL and -1 in NLL:
                NLL[1] = NLL[-1]

        c1 = ROOT.TCanvas("c1", "", 8, 30, 800, 800)
        NLL.zero()
        g = NLL.TGraph()
        mg.Add(g)

        g.SetLineColor(scan.color)
        g.SetLineStyle(scan.style)
        g.SetLineWidth(2)

        l.AddEntry(g, scan.title, "l")

    mg.Draw("AL")
    mg.GetXaxis().SetTitle(xtitle)
    mg.GetXaxis().SetRangeUser(xmin, xmax)
    mg.GetYaxis().SetTitle(ytitle)
    l.Draw()

    style.applycanvasstyle(c1)
    style.applyaxesstyle(mg)
    style.CMS(CMStext, luminosity, drawCMS=drawCMS)

    drawlines(CLtextposition, xmin=xmin, xmax=xmax)
    for ext in "png eps root pdf".split():
        outputfile = outputfile.replace("."+ext, "")
    for ext in "png eps root pdf".split():
        c1.SaveAs("{}.{}".format(outputfile, ext))

#https://root.cern.ch/phpBB3/viewtopic.php?t=10159
def GetNDC(x, y, logscale=False):
    ROOT.gPad.Update()#this is necessary!
    if logscale:
        y = math.log10(y)

    xndc = (x - ROOT.gPad.GetX1())/(ROOT.gPad.GetX2()-ROOT.gPad.GetX1())
    yndc = (y - ROOT.gPad.GetY1()) / (ROOT.gPad.GetY2() - ROOT.gPad.GetY1())
    return xndc, yndc

class XPos(object):
    def __init__(self, value):
        try:
            self.value = float(value)
            self.iscustom = False
            if not -1 <= self.value <= 1: raise ValueError
        except (TypeError, ValueError):
            if value in ("left", "right"):
                self.value = value
                self.iscustom = True
            else:
                raise ValueError("XPos argument has to be left, right, or a number from -1 to 1")
    def __eq__(self, other):
        return self.value == other
    def __ne__(self, other):
        return not self == other

    def TPaveText(self, ypos, logscale=False, yshift=.01, xsize=.1, ysize=.03):
        xsize = .1
        ysize = .03
        if self == "right":
            x2, y1 = GetNDC(.95, ypos, logscale)
            x1 = x2 - xsize
        elif self == "left":
            x1, y1 = GetNDC(-.95, ypos, logscale)
            x2 = x1 + xsize
        else:
            x1, y1 = GetNDC(self.value, ypos, logscale)
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
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--scan", nargs=4, action="append", help="example: --scan file1.root,file2.root,file3.root Expected kGreen+3 kDashed")
    parser.add_argument("outputfile", required=True)
    parser.add_argument("--legendposition", nargs=4, type=float, default=[.2, .7, .6, .9], help="--legendposition xmin xmax ymin ymax")
    parser.add_argument("--CLtextposition", default="left")
    parser.add_argument("--luminosity", required=True)
    parser.add_argument("--scanrange", nargs=2, type=float, default=[-1, 1], help="x axis limits")
    parser.add_argument("--POI", default="CMS_zz4l_fai1")
    parser.add_argument("--CMStext", default="Preliminary")
    parser.add_argument("--noCMS", dest="drawCMS", action="store_false", help="don't write CMS in the corner")
    parser.add_argument("--xtitle", required=True, help="x axis label")
    parser.add_argument("--ytitle", default="#minus2 #Deltaln L", help="y axis label")
    args = parser.parse_args()
    kwargs = args.__dict__

    scans = [Scan(filenames.split(","), title, parsecolor(color), parsestyle(style))
               for filenames, title, color, style in kwargs.pop("scan")]

    plotlimits(*scans, **kwargs)

    outputfile = kwargs["outputfile"]

    for ext in "png eps root pdf".split():
        outputfile = outputfile.replace("."+ext, "")
    with open(outputfile+".txt", 'w') as f:
        f.write("cd {} &&\n".format(os.getcwd()))
        f.write("python " + " ".join(sys.argv)+"\n")

        try:
            subprocess.check_call(["git", "status"])
        except subprocess.CalledProcessError:
            pass
        else:
            f.write("\n\n\n\n\n\ngit info:\n\n")
            f.write(subprocess.check_output(["git", "rev-parse", "HEAD"]))
            f.write("\n")
            f.write(subprocess.check_output(["git", "status"]))
            f.write("\n")
            f.write(subprocess.check_output(["git", "diff"]))
