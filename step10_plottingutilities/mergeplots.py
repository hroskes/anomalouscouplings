#!/usr/bin/env python
from glob import glob
import math
import os
import pipes
import subprocess
import sys

import ROOT

from Alignment.OfflineValidation.TkAlAllInOneTool.helperFunctions import replaceByMap

from helperstuff import config
from helperstuff.enums import Analysis, Production
from helperstuff.plotlimits import arrowatminimum, drawlines, xaxisrange
import helperstuff.stylefunctions as style
from helperstuff.utilities import cache, tfiles

from limits import findwhereyequals, Point

class Folder(object):
    def __init__(self, folder, title, color, analysis, subdir, plotname, graphnumber=None, repmap=None, linestyle=None, linewidth=None, secondcolumn=None):
        self.__folder, self.__title, self.color, self.analysis, self.subdir, self.plotname, self.graphnumber, self.linestyle, self.linewidth, self.secondcolumn = folder, title, color, Analysis(analysis), subdir, plotname, graphnumber, linestyle, linewidth, secondcolumn
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
        f = tfiles[replaceByMap(os.path.join(self.folder, self.plotname), self.repmap)]
        c = f.c1
        mg = c.GetListOfPrimitives()[1]
        graphs = mg.GetListOfGraphs()
        graphnumber = self.graphnumber
        if self.graphnumber is None:
            assert len(graphs) == 1
            graphnumber = 0
        graphs[graphnumber].SetLineColor(self.color)
        if self.linestyle is not None:
            graphs[graphnumber].SetLineStyle(self.linestyle)
        if self.linewidth is not None:
            graphs[graphnumber].SetLineWidth(self.linewidth)
        self.__xtitle = mg.GetXaxis().GetTitle()
        self.__ytitle = mg.GetYaxis().GetTitle()
        return graphs[graphnumber]
    @property
    def xtitle(self):
        self.graph
        return self.__xtitle
    @property
    def ytitle(self):
        self.graph
        return self.__ytitle
    def addtolegend(self, legend):
        legend.AddEntry(self.graph, self.title, "l")
        if self.secondcolumn is not None:
            legend.AddEntry(0, self.secondcolumn, "")

def mergeplots(analysis, **kwargs):
    drawlineskwargs = {}
    logscale = False
    PRL = False
    legendposition = (.2, .7, .6, .9)
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
        else:
            drawlineskwargs[kw] = kwarg

    analysis = Analysis(analysis)
    repmap = {"analysis": str(analysis)}
    subdir = ""
    plotname = "limit_lumi35.8671.root"
    folders = [
               Folder(".oO[analysis]Oo._allsysts", "Observed", 4, analysis, subdir, plotname="limit_lumi35.8671_7813_100,-1.0,1.0_100,-0.02,0.02.root", graphnumber=0, repmap=repmap, linestyle=1, linewidth=2),
               Folder(".oO[analysis]Oo._allsysts", "Expected", 4, analysis, subdir, plotname="limit_lumi35.8671_7813_100,-1.0,1.0_100,-0.02,0.02.root", graphnumber=1, repmap=repmap, linestyle=7, linewidth=2),
               Folder(".oO[analysis]Oo._allsysts", "Observed, 13 TeV", 1, analysis, subdir, plotname="limit_lumi35.8671_13_100,-1.0,1.0_100,-0.02,0.02.root", graphnumber=0, repmap=repmap, linestyle=1, linewidth=1),
               Folder(".oO[analysis]Oo._allsysts", "Expected, 13 TeV", 1, analysis, subdir, plotname="limit_lumi35.8671_13_100,-1.0,1.0_100,-0.02,0.02.root", graphnumber=1, repmap=repmap, linestyle=7, linewidth=1),
              ]
    outdir = ".oO[analysis]Oo._allsysts"

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

    c = ROOT.TCanvas("c1", "", 8, 30, 800, 800)
    mg.Draw("al")
    mg.GetXaxis().SetTitle(folders[0].xtitle)
    mg.GetYaxis().SetTitle(folders[0].ytitle)
    mg.GetXaxis().SetRangeUser(-1, 1)

    if PRL:
        mg.GetXaxis().CenterTitle()
        mg.GetYaxis().CenterTitle()

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
    style.CMS("", lumi=None, lumitext="5.1 fb^{{-1}} (7 TeV) + 19.7 fb^{{-1}} (8 TeV) + {:.1f} fb^{{-1}} (13 TeV)"
                                            .format(config.productionforcombine.dataluminosity+config.lumi2015))
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

if __name__ == "__main__":
    args = []
    kwargs = {}
    for arg in sys.argv[1:]:
        if "=" in arg:
            kwargs[arg.split("=")[0]] = arg.split("=", 1)[1]
        else:
            args.append(arg)
    if args or kwargs:
        mergeplots(*args, **kwargs)
    else:
        raise TypeError("Need to give args or kwargs")
