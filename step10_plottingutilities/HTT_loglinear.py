#!/usr/bin/env python
from array import array
from glob import glob
from itertools import izip
import math
import os
import pipes
import random
import shutil
import subprocess
import sys

import ROOT

from Alignment.OfflineValidation.TkAlAllInOneTool.helperFunctions import replaceByMap

from helperstuff import config
from helperstuff.enums import Analysis, Production
from helperstuff.plotlimits import arrowatminimum, drawlines, xaxisrange
import helperstuff.stylefunctions as style
from helperstuff.utilities import cache, cd, mkdtemp, PlotCopier, tfiles

from mergeplots import Folder

analyses = "fa3", "fa2", "fL1", "fL1Zg"
setmax = 1
def saveasdir(forWIN=False):
    return os.path.join(config.plotsbasedir, "limits", "HTT")
baseplotname = "limit_lumi80.15.root"

def applystyle(mgs, mglogs, folders, xdivides, ydivide):
    assert len(mgs) == len(mglogs) == len(xdivides)+1
    for mg, mglog, xmin, xmax in izip(mgs, mglogs, [-setmax]+list(xdivides), list(xdivides)+[setmax]):
        mglog.GetXaxis().SetTitle(folders[0].xtitle)
        mglog.GetXaxis().SetRangeUser(xmin, xmax)
        mglog.GetXaxis().CenterTitle()
        mglog.GetYaxis().CenterTitle()
        mglog.SetMinimum(ydivide)
        mglog.SetMaximum(250)
        style.applyaxesstyle(mglog)
        mglog.GetXaxis().SetLabelOffset(9999999)
        mglog.GetXaxis().SetTitleOffset(9999999)
        mglog.GetYaxis().SetLabelSize(.1)
        mglog.GetYaxis().SetTitleSize(.1)

        mg.GetXaxis().SetTitle(folders[0].xtitle)
        mg.GetXaxis().SetRangeUser(xmin, xmax)
        mg.GetXaxis().CenterTitle()
        mg.GetYaxis().CenterTitle()
        mg.SetMinimum(0)
        mg.SetMaximum(ydivide)
        style.applyaxesstyle(mg)
        mg.GetXaxis().SetLabelSize(.1)
        mg.GetYaxis().SetLabelSize(.1)
        mg.GetXaxis().SetTitleSize(.12)
        mg.GetYaxis().SetTitleSize(.1)

def getplotname(analysis):
    return "{}_{}".format(analysis, baseplotname)


def PRL_loglinear(**kwargs):
    commondrawlineskwargs = {
                             "logscale": False,  #the lines are in the linear part
                             "xsize": .2,
                             "ysize": .045,
                             "textsize": .09,
                             "yshift68": .03,
                             "yshift95": .03,
                            }
    markerposition = None
    onlyanalysis = None
    ydivide = 11
    xdivides = -.05, .05
    saveas = None
    forWIN = False
    for kw, kwarg in kwargs.iteritems():
        if kw == "ydivide":
            ydivide = float(kwarg)
        elif kw == "xdivides":
            xdivides = sorted(float(_) for _ in kwarg)
            assert len(xdivides) == 2, xdivides
        elif kw == "markerposition":
            markerposition = kwarg
        elif kw == "analysis":
            onlyanalysis = kwarg
        elif kw == "saveas":
            saveas = kwarg
        elif kw == "forWIN":
            forWIN = kwarg
        else:
            commondrawlineskwargs[kw] = kwarg

    for k, v in commondrawlineskwargs.items():
        if k == "xpostext":
            try:
                commondrawlineskwargs[k] = float(v)
            except ValueError:
                pass
        elif k in ("xmin", "xmax"):
            commondrawlineskwargs[k] = float(v)

    for i, (analysis, letter) in reversed(list(enumerate(zip(analyses, "abcd"), start=1))):
        if analysis != onlyanalysis is not None: continue
        if analysis == "fa3":
            x1legend = .32
            CLtextposition=-.9
        elif analysis == "fa2":
            x1legend = .5
            CLtextposition=-.9
        elif analysis == "fL1":
            x1legend = .15
            CLtextposition=.65
        elif analysis == "fL1Zg":
            x1legend = .15
            CLtextposition=-.9
        else:
            assert False
        legendposition = x1legend, .3, x1legend+.45, .8

        c = plotcopier.TCanvas("c{}".format(random.randint(1, 1000000)), "", 8, 30, 1600, 1600)

        leftmargin = .2
        rightmargin = .045 #apply to the individual pads or 1 of the x axis gets cut off
        topmargin = .02
        bottommargin = .225
        assert abs((leftmargin + rightmargin) - (topmargin + bottommargin)) < 1e-5, (leftmargin + rightmargin, topmargin + bottommargin)
        c.SetLeftMargin(leftmargin)
        c.SetRightMargin(0)
        c.Divide(3, 2, 0, 0)
        logpads = [c.cd(1), c.cd(2), c.cd(3)]
        linearpads = [c.cd(4), c.cd(5), c.cd(6)]
        legendpad = logpads[1]
        for _ in linearpads + logpads:
            _.SetTicks()

        analysis = Analysis(analysis)
        repmap = {"analysis": str(analysis)}
        subdir = ""
        folders = [
                   Folder(".oO[analysis]Oo._August17combination", "Observed", 2, analysis, subdir, plotname="limit_lumi77.45_7813_101,-1.0,1.0_101,-0.02,0.02.root", graphnumber=0, repmap=repmap, linestyle=1, linewidth=2),
                   Folder(".oO[analysis]Oo._August17combination", "Expected", 2, analysis, subdir, plotname="limit_lumi77.45_7813_101,-1.0,1.0_101,-0.02,0.02.root", graphnumber=-1, repmap=repmap, linestyle=7, linewidth=2),
                   Folder(".oO[analysis]Oo._August17combination", "Observed, 2016+2017", 1, analysis, subdir, plotname="limit_lumi77.45_101,-1.0,1.0_101,-0.02,0.02.root", graphnumber=0, repmap=repmap, linestyle=1, linewidth=1),
                   Folder(".oO[analysis]Oo._August17combination", "Expected, 2016+2017", 1, analysis, subdir, plotname="limit_lumi77.45_101,-1.0,1.0_101,-0.02,0.02.root", graphnumber=-1, repmap=repmap, linestyle=7, linewidth=1),
                  ]

        mg = ROOT.TMultiGraph("limit", "")
        for folder in folders:
            mg.Add(folder.graph)

        drawlineskwargs = commondrawlineskwargs.copy()
        drawlineskwargs["xpostext"] = CLtextposition
        drawlineskwargs["arbitraryparameter"] = analysis, markerposition

        if markerposition:
            markerposition = [array('d', [_]) for _ in markerposition]
            marker = ROOT.TGraph(1, *markerposition)
            marker.SetMarkerStyle(20)
            marker.SetMarkerColor(4)
            marker.SetMarkerSize(3)
            if marker.GetY()[0] > 0:
                mg.Add(marker, "P")

        mglog = mg.Clone()
        if markerposition and marker.GetY()[0] <= 0:
            mg.Add(marker, "P")
        mgs = [mg, mg.Clone(), mg.Clone()]
        mglogs = [mglog, mglog.Clone(), mglog.Clone()]
        for logpad, mglog in izip(logpads, mglogs):
            logpad.cd()
            logpad.SetLogy()
            mglog.Draw("al")
            logpad.SetTopMargin(topmargin*2)
        logpads[-1].SetRightMargin(rightmargin)
        style.subfig(letter, textsize=.11, x1=.87, x2=.91, y1=.87, y2=.91)

        for linearpad, mg in izip(linearpads, mgs):
            linearpad.cd()
            mg.Draw("al")
            linearpad.SetBottomMargin(bottommargin*2)
        linearpads[-1].SetRightMargin(rightmargin)

        applystyle(mgs, mglogs, folders, xdivides, ydivide)

        drawlines(**drawlineskwargs)

        legendpad.cd()
        l = ROOT.TLegend(*legendposition)
        l.SetBorderSize(0)
        l.SetFillStyle(0)
        for folder in folders:
            folder.addtolegend(l)
        if all(folder.secondcolumn is not None for folder in folders):
            l.SetNColumns(2)
        l.Draw()
        c.cd()
        style.applycanvasstyle(c)
        style.CMS("", lumi=None, lumitext="5.1 fb^{{-1}} (7 TeV) + 19.7 fb^{{-1}} (8 TeV) + {:.1f} fb^{{-1}} (13 TeV)"
                                                .format(sum(_.dataluminosity for _ in config.productionsforcombine)+config.lumi2015),
                      x1=0.007, x2=1.01, #???
                      drawCMS=False, extratextsize=.039)
        if forWIN:
            if analysis == "fa2":
                style.CMS("", x1=0.12, x2=1.025, y1=.86, y2=.94, CMStextsize=.06, extratextsize=.039)
                style.CMS("Preliminary", x1=0.0, x2=1.025, y1=.8, y2=.86, CMStextsize=.06, extratextsize=.039, drawCMS=False)
            else:
                style.CMS("", x1=0.12, x2=1.025, y1=.86, y2=.94, CMStextsize=.06, extratextsize=.039)
                style.CMS("Preliminary", x1=0.15, x2=1.025, y1=.85, y2=.93, CMStextsize=.06, extratextsize=.039, drawCMS=False)
        else:
            style.CMS("", x1=0.12, x2=1.025, y1=.86, y2=.94, CMStextsize=.06)
        yaxislabel(folders[0].ytitle).Draw()

        try:
            os.makedirs(saveasdir(forWIN))
        except OSError:
            pass
        plotname = getplotname(analysis)
        assert ".root" in plotname
        if saveas is not None:
            c.SaveAs(saveas)
        else:
            for ext in "png eps root pdf".split():
                c.SaveAs(os.path.join(saveasdir(forWIN), replaceByMap(plotname.replace("root", ext), repmap)))
            with plotcopier.open(os.path.join(saveasdir(forWIN), replaceByMap(plotname.replace("root", "txt"), repmap)), "w") as f:
                f.write(" ".join(["python"]+[pipes.quote(_) for _ in sys.argv]))
                f.write("\n\n\n\n\n\ngit info:\n\n")
                f.write(subprocess.check_output(["git", "rev-parse", "HEAD"]))
                f.write("\n")
                f.write(subprocess.check_output(["git", "status"]))
                f.write("\n")
                f.write(subprocess.check_output(["git", "diff"]))
            if hasattr(config, "svndir"):
                shutil.copyfile(
                                os.path.join(saveasdir(forWIN), replaceByMap(plotname.replace("root", "pdf"), repmap)),
                                os.path.join(config.svndir, "papers", "HIG-17-011", "trunk", "Figures", "fig3{}.pdf".format(letter))
                               )

def animations(**kwargs):
    from projections import Projections
    forWIN = kwargs.get("forWIN", False)
    for analysis in analyses:
        tmpdir = mkdtemp()
        convertcommand = ["gm", "convert", "-loop", "0"]
        animation = Projections.animationstepsforniceplots(analysis)
        lastdelay = None
        for i, step in enumerate(animation):
            if step.delay != lastdelay:
                convertcommand += ["-delay", str(step.delay)]
            convertcommand += ["-trim", os.path.join(tmpdir, "{}.pdf".format(i))]
            PRL_loglinear(
                          analysis=analysis,
                          saveas=os.path.join(tmpdir, "{}.pdf".format(i)),
                          markerposition=(step.fai_decay, step.deltaNLL),
                          **kwargs
                         )

        finalplot = os.path.join(saveasdir(forWIN), getplotname(analysis).replace("root", "gif"))
        convertcommand.append(finalplot)
        #http://stackoverflow.com/a/38792806/5228524
        #subprocess.check_call(convertcommand)
        os.system(" ".join(pipes.quote(_) for _ in convertcommand))
        shutil.rmtree(tmpdir)


@cache
def yaxislabel(label, textsize=.06):
    pt = ROOT.TPaveText(.01, 0, .03, 1, "brNDC")
    pt.SetBorderSize(0)
    pt.SetFillStyle(0)
    pt.SetTextAlign(22)
    pt.SetTextFont(42)
    pt.SetTextSize(.2)
    text = pt.AddText(.5,.5,label)
    text.SetTextSize(textsize)
    text.SetTextAngle(90)
    return pt

if __name__ == "__main__":
    args = []
    kwargs = {}
    for arg in sys.argv[1:]:
        if "=" in arg:
            kwargs[arg.split("=")[0]] = arg.split("=", 1)[1]
        else:
            args.append(arg)
    function = PRL_loglinear
    if args and args[0] == "animations":
        args = args[1:]
        function = animations
    with PlotCopier() as plotcopier:
        function(*args, **kwargs)
