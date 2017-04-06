#!/usr/bin/env python
from glob import glob
import os
import pipes
import subprocess
import sys

import ROOT

from Alignment.OfflineValidation.TkAlAllInOneTool.helperFunctions import replaceByMap

from helperstuff import config
from helperstuff.enums import Analysis, Production
from helperstuff.plotlimits import drawlines, xaxisrange
import helperstuff.stylefunctions as style
from helperstuff.utilities import cache, tfiles

class Folder(object):
    def __init__(self, folder, title, color, analysis, subdir, plotname, graphnumber=None, repmap=None):
        self.__folder, self.__title, self.color, self.analysis, self.subdir, self.plotname, self.graphnumber = folder, title, color, Analysis(analysis), subdir, plotname, graphnumber
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

def mergeplots(analysis, **kwargs):
    drawlineskwargs = {}
    logscale = False
    for kw, kwarg in kwargs.iteritems():
        if kw == "logscale":
            logscale = bool(int(kwarg))
            drawlineskwargs[kw] = kwarg
        else:
            drawlineskwargs[kw] = kwarg

    analysis = Analysis(analysis)
    repmap = {"analysis": str(analysis)}
    subdir = ""
    plotname = "limit_lumi35.8671.root"
    folders = [
               Folder(".oO[analysis]Oo._allsysts", "Observed", 1, analysis, subdir, plotname="limit_lumi35.8671_7813_100,-1.0,1.0_100,-0.02,0.02.root", graphnumber=0, repmap=repmap),
               Folder(".oO[analysis]Oo._allsysts", "Expected", 1, analysis, subdir, plotname="limit_lumi35.8671_7813_100,-1.0,1.0_100,-0.02,0.02.root", graphnumber=1, repmap=repmap),
               Folder(".oO[analysis]Oo._allsysts", "Observed, 13 TeV", 6, analysis, subdir, plotname="limit_lumi35.8671_13_100,-1.0,1.0_100,-0.02,0.02.root", graphnumber=0, repmap=repmap),
               Folder(".oO[analysis]Oo._allsysts", "Expected, 13 TeV", 6, analysis, subdir, plotname="limit_lumi35.8671_13_100,-1.0,1.0_100,-0.02,0.02.root", graphnumber=1, repmap=repmap),
              ]
    outdir = ".oO[analysis]Oo._allsysts"

    mg = ROOT.TMultiGraph("limit", "")
    if analysis == "fa3":
        l = ROOT.TLegend(.3, .6, .6, .9)
    else:
        l = ROOT.TLegend(.6, .6, .9, .9)
    l.SetBorderSize(0)
    l.SetFillStyle(0)
    for folder in folders:
        mg.Add(folder.graph)
        folder.addtolegend(l)

    c = ROOT.TCanvas("c", "", 8, 30, 800, 800)
    mg.Draw("al")
    mg.GetXaxis().SetTitle(folders[0].xtitle)
    mg.GetYaxis().SetTitle(folders[0].ytitle)
    mg.GetXaxis().SetRangeUser(-1, 1)

    if logscale:
        c.SetLogy()
        mg.SetMinimum(0.1)
        plotname = plotname.replace(".root", "_log.root")

    l.Draw()
    style.applycanvasstyle(c)
    style.applyaxesstyle(mg)
    style.CMS("Preliminary", Production(config.productionsforcombine[0]).dataluminosity)
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
