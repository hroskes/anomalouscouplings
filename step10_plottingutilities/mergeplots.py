#!/usr/bin/env python
if __name__ == "__main__":
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
    def __init__(self, folder, title, color, analysis, subdir, plotname, graphnumber=None, repmap=None, linestyle=None, linewidth=None, secondcolumn=None, removepoints=None, forcepoints=None):
        if removepoints is None: removepoints = []
        if forcepoints is None: forcepoints = {}
        self.__folder, self.__title, self.color, self.analysis, self.subdir, self.plotname, self.graphnumber, self.linestyle, self.linewidth, self.secondcolumn, self.removepoints, self.forcepoints = folder, title, color, Analysis(analysis), subdir, plotname, graphnumber, linestyle, linewidth, secondcolumn, removepoints, forcepoints
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
            for forcex, forcey in self.forcepoints.iteritems():
                forceindices = np.isclose(x, forcex)
                assert sum(forceindices) == 1, forceindices
                print y[forceindices], forcey
                y[forceindices] = forcey
                print y[forceindices]
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
    zoom = kwargs.pop("zoom", "zoom")
    legendposition = kwargs.pop("legendposition", {
      "zoom": {
        Analysis("fa3"): (.2, .5, .45, .65),
        Analysis("fa2"): (.42, .75, .67, .9),
        Analysis("fL1"): (.2, .75, .45, .9),
        Analysis("fL1Zg"): (.4, .75, .65, .9),
      },
      "": {
        Analysis("fa3"): (.3, .75, .55, .9),
        Analysis("fa2"): (.2, .75, .45, .9),
        Analysis("fL1"): (.2, .75, .45, .9),
        Analysis("fL1Zg"): (.2, .75, .45, .9),
      },
      "medium": {
        Analysis("fa3"): (.2, .75, .45, .9),
        Analysis("fa2"): (.2, .6, .45, .75),
        Analysis("fL1"): (.2, .75, .45, .9),
        Analysis("fL1Zg"): (.2, .75, .45, .9),
      },
    }[zoom][analysis])
    ymax = kwargs.pop("ymax", 40 if zoom else 1000)
    defaultxmax = {"zoom": 0.0055, "": 1, "medium": 0.1}[zoom]
    xmin, xmax = kwargs.pop("xmin", -defaultxmax), kwargs.pop("xmax", defaultxmax)
    drawlineskwargs = {"xpostext": kwargs.pop("CLtextposition", {"zoom": 0.003, "": 99, "medium": 0.03 if analysis != "fL1Zg" else -0.09}[zoom]), "xmin": xmin, "xmax": xmax}
    subdir = kwargs.pop("subdir", "")
    lumi = kwargs.pop("lumi", None)
    outdir = kwargs.pop("outdir", "fa3fa2fL1fL1Zg_morecategories_writeup")
    plotcopier = kwargs.pop("plotcopier", ROOT)
    drawlineskwargs.update(kwargs)

    repmap = {"analysis": str(analysis)}
    plotname = "limit_lumi3000.00_scan.oO[analysis]Oo._compare"+("_"+zoom if zoom else "")+".root"
    if analysis == "fa3":
        removepoints1 = [-0.6, 0.6]
        removepoints2 = []
        removepoints3 = []
    if analysis == "fa2":
        removepoints1 = [.22, .24, .3]
        removepoints2 = []
        removepoints3 = []
    if analysis == "fL1":
        removepoints1 = [-.32, -.3, -.18, .3, -.52, -.5]
        removepoints2 = [-.1, -.06, .06, .08, .24, .66]
        removepoints3 = []
    if analysis == "fL1Zg":
        removepoints1 = [-.68, -.66, -.64, -.6, -0.56, -0.54, -0.52]
        removepoints2 = []
        removepoints3 = []
    folders = [
               Folder("fa3fa2fL1fL1Zg_morecategories_writeup/", "MELA", 2, analysis, subdir, plotname="limit_lumi3000.00_scan.oO[analysis]Oo._101,-0.02,0.02_merged.root", graphnumber=0, repmap=repmap, linestyle=2, linewidth=2, removepoints=removepoints1),
               Folder("fa3fa2fL1fL1Zg_STXS_writeup/", "STXS stage 1", 4, analysis, subdir, plotname="limit_lumi3000.00_scan.oO[analysis]Oo._101,-0.02,0.02_merged.root", graphnumber=0, repmap=repmap, linestyle=2, linewidth=2, removepoints=removepoints2),
               Folder("fa3fa2fL1fL1Zg_decay_writeup/", "decay only", ROOT.kGreen+3, analysis, subdir, plotname="limit_lumi3000.00_scan.oO[analysis]Oo._merged.root", graphnumber=0, repmap=repmap, linestyle=2, linewidth=2, removepoints=removepoints3),
              ]

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
    mg.GetXaxis().SetLimits(xmin, xmax)
    mg.SetMinimum(0)

    if ymax is not None:
        mg.SetMaximum(ymax)

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
    if zoom != "": drawlines(**drawlineskwargs)
    saveasdir = replaceByMap(os.path.join(config.plotsbasedir, "limits", subdir, outdir), repmap)
    try:
        os.makedirs(saveasdir)
    except OSError:
        pass
    assert ".root" in plotname
    for ext in "png eps root pdf C".split():
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
            for zoom in "zoom", "", "medium":
                mergeplots(analysis, plotcopier=pc, zoom=zoom, **args.__dict__)
