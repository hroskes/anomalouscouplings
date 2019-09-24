#!/usr/bin/env python
import argparse
p = argparse.ArgumentParser()
args = p.parse_args()

from glob import glob
import math
import os
import pipes
import subprocess
import sys

import numpy as np

import ROOT

from Alignment.OfflineValidation.TkAlAllInOneTool.helperFunctions import replaceByMap

from helperstuff import config
from helperstuff.enums import Analysis, Production
from helperstuff.plotlimits import arrowatminimum, drawlines, xaxisrange
import helperstuff.stylefunctions as style
from helperstuff.utilities import cache, PlotCopier, TFile

from limits import findwhereyequals, Point

class Folder(object):
    def __init__(self, folder, title, color, analysis, subdir, plotname, graphnumber=None, repmap=None, linestyle=None, linewidth=None, secondcolumn=None, removepoints=None):
        if removepoints is None: removepoints = []
        self.__folder, self.__title, self.color, self.analysis, self.subdir, self.plotname, self.graphnumber, self.linestyle, self.linewidth, self.secondcolumn, self.removepoints = folder, title, color, Analysis(analysis), subdir, plotname, graphnumber, linestyle, linewidth, secondcolumn, removepoints
        self.repmap = {
                       "analysis": str(self.analysis),
                      }
        if repmap is not None: self.repmap.update(repmap)
    @property
    def folder(self):
        result = os.path.join(config.plotsbasedir, "limits", self.subdir, replaceByMap(self.__folder, self.repmap))
        gl = glob(result)
        if not gl:
            raise ValueError("{} does not exist!".format(result))
        if len(gl) > 1:
            raise ValueError("{} returns more than one match!".format(result))
        return gl[0]
    @property
    def title(self):
        return replaceByMap(self.__title, self.repmap)
    @title.setter
    def title(self, newtitle):
        self.__title = newtitle
    @property
    @cache
    def graph(self):
        with TFile(replaceByMap(os.path.join(self.folder, self.plotname), self.repmap)) as f:
            c = f.c1
            mg = c.GetListOfPrimitives()[1]
            graphs = mg.GetListOfGraphs()
            graphnumber = self.graphnumber
            if self.graphnumber is None:
                assert len(graphs) == 1
                graphnumber = 0
            self.__xtitle = mg.GetXaxis().GetTitle()
            self.__ytitle = mg.GetYaxis().GetTitle()
            g = graphs[graphnumber]
            n = g.GetN()
            x = np.array([g.GetX()[i] for i in xrange(n) if not any(np.isclose(g.GetX()[i], _) for _ in self.removepoints)])
            y = np.array([g.GetY()[i] for i in xrange(n) if not any(np.isclose(g.GetX()[i], _) for _ in self.removepoints)])
            print self.removepoints, g.GetN(), len(x)
            y -= min(y)
            newg = ROOT.TGraph(len(x), x, y)
            newg.SetLineColor(self.color)
            if self.linestyle is not None:
                newg.SetLineStyle(self.linestyle)
            if self.linewidth is not None:
                newg.SetLineWidth(self.linewidth)
            return newg
    @property
    def xtitle(self):
        self.graph
        return self.__xtitle
    @property
    def ytitle(self):
        self.graph
        return self.__ytitle
    def addtolegend(self, legend, graph=None):
        if graph is None: graph = self.graph
        legend.AddEntry(self.graph, self.title, "l")
        if self.secondcolumn is not None:
            legend.AddEntry(0, self.secondcolumn, "")

def mergeplots(analysis, **kwargs):
    analysis = Analysis(analysis)
    drawlineskwargs = {"xpostext": 0.0025, "xmin": -0.005, "xmax": 0.005}
    logscale = False
    PRL = False
    legendposition = {
      Analysis("fa3"): (.2, .5, .45, .65),
      Analysis("fa2"): (.42, .75, .67, .9),
      Analysis("fL1"): (.2, .75, .45, .9),
      Analysis("fL1Zg"): (.4, .75, .65, .9),
    }[analysis]
    ymax = 40
    xmin, xmax = -0.005, 0.005
    subdir = ""
    lumi = None
    outdir = "fa3fa2fL1fL1Zg_morecategories_writeup"
    plotcopier = ROOT
    for kw, kwarg in kwargs.iteritems():
        if kw == "logscale":
            logscale = bool(int(kwarg))
            drawlineskwargs[kw] = kwarg
        elif kw == "PRL":
            PRL = bool(int(kwarg))
            drawlineskwargs[kw] = kwarg
        elif kw == "CLtextposition":
            drawlineskwargs["xpostext"] = kwarg
        elif kw == "legendposition":
            try:
                legendposition = [float(a) for a in kwarg.split(",")]
                if len(legendposition) != 4: raise ValueError
            except ValueError:
                raise ValueError("legendposition has to contain 4 floats separated by commas!")
        elif kw == "ymax":
            ymax = float(kwarg)
        elif kw == "xmin":
            xmin = float(kwarg)
            drawlineskwargs[kw] = kwarg
        elif kw == "xmax":
            xmax = float(kwarg)
            drawlineskwargs[kw] = kwarg
        elif kw == "subdir":
            subdir = kwarg
        elif kw == "outdir":
            outdir = kwarg
        elif kw == "lumi":
            lumi = float(kwarg)
        elif kw == "plotcopier":
            plotcopier = kwarg
        else:
            drawlineskwargs[kw] = kwarg

    repmap = {"analysis": str(analysis)}
    plotname = "limit_lumi3000.00_scan.oO[analysis]Oo._compare.root"
    folders = [
               Folder("fa3fa2fL1fL1Zg_morecategories_writeup/", "MELA", 2, analysis, subdir, plotname="limit_lumi3000.00_scan.oO[analysis]Oo._101,-0.02,0.02_merged.root", graphnumber=0, repmap=repmap, linestyle=2, linewidth=2),
               Folder("fa3fa2fL1fL1Zg_STXS_writeup/", "STXS stage 1", 4, analysis, subdir, plotname="limit_lumi3000.00_scan.oO[analysis]Oo._101,-0.02,0.02_merged.root", graphnumber=0, repmap=repmap, linestyle=2, linewidth=2),
              ]

    if logscale and config.minimainlegend:
        for folder in folders:
            graph = folder.graph
            folder.title = folder.title.replace("Observed, 13 TeV", "Obs. 13 TeV")
            folder.title = folder.title.replace("Expected, 13 TeV", "Exp. 13 TeV")
            minimum = Point(float("nan"), float("inf"))
            CLleft = CLright = None
            for i, x, y in zip(xrange(graph.GetN()), graph.GetX(), graph.GetY()):
                if y > 1:
                    if CLright is None and i > 0:
                        CLright = findwhereyequals(1, Point(lastx, lasty), Point(x, y))
                        print "right", lastx, "    ", lasty, "    ", x, "    ", y
                if y < 1:
                    CLright = None
                    if CLleft is None:
                        CLleft = findwhereyequals(1, Point(lastx, lasty), Point(x, y))
                        print "left ", lastx, "    ", lasty, "    ", x, "    ", y
                if y < minimum.y:
                    minimum = Point(x, y)

                lasti, lastx, lasty = i, x, y

            minimum = minimum.x
            CLplus = CLright - minimum
            CLminus = minimum - CLleft
            if analysis == "fa3" and "Exp" in folder.title: CLplus = CLminus = (CLplus+CLminus)/2
            ndigits = min(1 - math.floor(math.log10(max(CLplus, CLminus))), 3)
            appendfmt = "{:.%df}^{{{:+.%df}}}_{{{:+.%df}}}" % (ndigits, ndigits, ndigits)
            if float(("{:.%df}"%ndigits).format(minimum)) == 0:
                minimum = 0 #to turn -0.00 --> 0.00
            append = appendfmt.format(minimum, CLplus, -CLminus)
            append = append.replace("+", "#plus").replace("-", "#minus")
            folder.secondcolumn = append
            print folder.secondcolumn


    mg = ROOT.TMultiGraph("limit", "")
    l = ROOT.TLegend(*legendposition)
    l.SetBorderSize(0)
    l.SetFillStyle(0)
    for folder in folders:
        mg.Add(folder.graph)
        folder.addtolegend(l)
    if all(folder.secondcolumn is not None for folder in folders):
        l.SetNColumns(2)

    c = plotcopier.TCanvas("c1", "", 8, 30, 800, 800)
    mg.Draw("al")
    mg.GetXaxis().SetTitle(folders[0].xtitle.replace("a2", "g2").replace("a3", "g4"))
    mg.GetYaxis().SetTitle(folders[0].ytitle)
    mg.GetXaxis().SetRangeUser(xmin, xmax)
    mg.SetMinimum(0)

    if PRL:
        mg.GetXaxis().CenterTitle()
        mg.GetYaxis().CenterTitle()

    if ymax is not None:
        if logscale: raise ValueError("can't set ymax and logscale!")
        mg.SetMaximum(ymax)

    if logscale:
        c.SetLogy()
        mg.SetMinimum(0.1)
        mg.SetMaximum(120)
        plotname = plotname.replace(".root", "_log.root")
        for folder in folders:
            if config.arrowsatminima:
                if "Observed" in folder.title:
                    if analysis == "fa3" and folder.title == "Observed": continue
                    if analysis in ("fa3", "fL1Zg"):
                        abovexaxis = False
                    else:
                        abovexaxis = True
                    arrowatminimum(folder.graph, abovexaxis=abovexaxis).Draw()

    l.Draw()
    style.applycanvasstyle(c)
    style.applyaxesstyle(mg)
    if lumi is not None:
        style.CMS("Preliminary", lumi=None, lumitext="{:.1f} fb^{{-1}} (13 TeV)".format(lumi))
    for k, v in drawlineskwargs.items():
        if k == "xpostext":
            try:
                drawlineskwargs[k] = float(v)
            except ValueError:
                pass
        elif k in ("xmin", "xmax"):
            drawlineskwargs[k] = float(v)
    drawlines(**drawlineskwargs)
    saveasdir = replaceByMap(os.path.join(config.plotsbasedir, "limits", subdir, outdir), repmap)
    try:
        os.makedirs(saveasdir)
    except OSError:
        pass
    assert ".root" in plotname
    for ext in "png eps root pdf".split():
        c.SaveAs(os.path.join(saveasdir, replaceByMap(plotname.replace("root", ext), repmap)))
    with open(os.path.join(saveasdir, replaceByMap(plotname.replace("root", "txt"), repmap)), "w") as f:
        f.write(" ".join(["python"]+[pipes.quote(_) for _ in sys.argv]))
        f.write("\n\n")
        if replaceByMap(outdir, repmap).startswith(str(analysis)+"_"):
            restofoutdir = replaceByMap(outdir, repmap).replace(str(analysis)+"_", "")
            f.write("python limits.py {} {} ".format(analysis, restofoutdir))
            if subdir: f.write("subdirectory={} ".format(pipes.quote(subdir)))
            f.write("plotname="+plotname+" ")
        f.write("\n\n\n\n\n\ngit info:\n\n")
        f.write(subprocess.check_output(["git", "rev-parse", "HEAD"]))
        f.write("\n")
        f.write(subprocess.check_output(["git", "status"]))
        f.write("\n")
        f.write(subprocess.check_output(["git", "diff"]))

    if plotcopier != ROOT:
        pc.copy(os.path.join(saveasdir, replaceByMap(plotname.replace("root", "txt"), repmap)))

if __name__ == "__main__":
    with PlotCopier() as pc:
        for analysis in "fa3", "fa2", "fL1", "fL1Zg":
            mergeplots(analysis, plotcopier=pc, **args.__dict__)
