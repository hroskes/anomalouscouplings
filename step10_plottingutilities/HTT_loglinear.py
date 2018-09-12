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
    return os.path.join(config.plotsbasedir, "limits")
baseplotname = "limit_lumi80.15.root"

def applystyle(mgs, mglogs, folders, xboundaries, xdivides, ydivide):
    assert len(mgs) == len(mglogs) == len(xdivides)+1
    for mg, mglog, xmin, xmax, fractionxmin, fractionxmax in izip(mgs, mglogs, [-setmax]+list(xdivides), list(xdivides)+[setmax], xboundaries[:-1], xboundaries[1:]):
        mglog.GetXaxis().SetLimits(xmin, xmax)
        mglog.GetXaxis().CenterTitle()
        mglog.GetYaxis().CenterTitle()
        mglog.SetMinimum(ydivide)
        mglog.SetMaximum(250)
        style.applyaxesstyle(mglog)
        mglog.GetXaxis().SetLabelOffset(9999999)
        mglog.GetXaxis().SetTitleOffset(9999999)
        mglog.GetYaxis().SetLabelSize(.1)
        mglog.GetYaxis().SetTitleSize(.1)

        mg.GetXaxis().SetLimits(xmin, xmax)
        mg.GetXaxis().CenterTitle()
        mg.GetYaxis().CenterTitle()
        mg.SetMinimum(0)
        mg.SetMaximum(ydivide)
        style.applyaxesstyle(mg)
        mg.GetXaxis().SetLabelSize(.1 / (3*(fractionxmax - fractionxmin)))
        mg.GetXaxis().SetTitleOffset(0.7)
        mg.GetYaxis().SetLabelSize(.1 / (3*(fractionxmax - fractionxmin)))
        mg.GetXaxis().SetTitleSize(.15 / (3*(fractionxmax - fractionxmin)))
        mg.GetYaxis().SetTitleSize(.1 / (3*(fractionxmax - fractionxmin)))
        mg.GetXaxis().SetLabelOffset(-0.0165)

    mgs[0].GetXaxis().SetLabelOffset(0.007)
    mgs[-1].GetXaxis().SetLabelOffset(-0.012)

    mgs[len(mgs)/2].GetXaxis().SetTitle(folders[0].xtitle)
    for mg, mglog in izip(mgs[1:], mglogs[1:]):
        mglog.GetYaxis().SetLabelOffset(9999999)
        mglog.GetYaxis().SetTitleOffset(9999999)
        mg.GetYaxis().SetLabelOffset(9999999)
        mg.GetYaxis().SetTitleOffset(9999999)

def getplotname(analysis):
    return "{}_{}".format(analysis, baseplotname)


def PRL_loglinear(**kwargs):
    commondrawlineskwargs = {
                             "logscale": False,  #the lines are in the linear part
                             "xsize": .2,
                             "ysize": .045,
                             "yshift68": .03,
                             "yshift95": .03,
                            }
    markerposition = None
    onlyanalysis = None
    ydivide = 11
    xdivides = -.03, .03
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
            CLtextposition=-.8
        elif analysis == "fa2":
            CLtextposition=-.8
        elif analysis == "fL1":
            CLtextposition=.65
        elif analysis == "fL1Zg":
            CLtextposition=-.8
        else:
            assert False

        c = plotcopier.TCanvas("c{}".format(random.randint(1, 1000000)), "", 8, 30, 1600, 1600)

        leftmargin = .1
        rightmargin = .02 #apply to the individual pads or 1 of the x axis gets cut off
        topmargin = .07
        bottommargin = .12
        #assert abs((leftmargin + rightmargin) - (topmargin + bottommargin)) < 1e-5, (leftmargin + rightmargin, topmargin + bottommargin)
        c.SetLeftMargin(0)
        c.SetRightMargin(0)
        c.SetTopMargin(0)
        c.SetBottomMargin(0)
        xboundaries = [0, leftmargin+(1-leftmargin-rightmargin)/3, leftmargin+(1-leftmargin-rightmargin)*2/3, 1]
        yboundaries = [0, bottommargin+(1-bottommargin-topmargin)/2, 1]
        legendposition = 0, .2, 1, .8*(1-topmargin)
        logpads = [
          ROOT.TPad(c.GetName()+"_1", "_1", xboundaries[0], yboundaries[1], xboundaries[1], yboundaries[2]),
          ROOT.TPad(c.GetName()+"_2", "_2", xboundaries[1], yboundaries[1], xboundaries[2], yboundaries[2]),
          ROOT.TPad(c.GetName()+"_3", "_3", xboundaries[2], yboundaries[1], xboundaries[3], yboundaries[2]),
        ]
        linearpads = [
          ROOT.TPad(c.GetName()+"_4", "_4", xboundaries[0], yboundaries[0], xboundaries[1], yboundaries[1]),
          ROOT.TPad(c.GetName()+"_5", "_5", xboundaries[1], yboundaries[0], xboundaries[2], yboundaries[1]),
          ROOT.TPad(c.GetName()+"_6", "_6", xboundaries[2], yboundaries[0], xboundaries[3], yboundaries[1]),
        ]
        legendpad = logpads[1]
        for _ in linearpads + logpads:
            _.Draw()
            _.SetTicks()
            _.SetLeftMargin(0)
            _.SetRightMargin(0)
            _.SetTopMargin(0)
            _.SetBottomMargin(0)
        logpads[-1].SetRightMargin(rightmargin*len(logpads))

        analysis = Analysis(analysis)
        repmap = {"analysis": str(analysis)}
        subdir = ""
        folders = [
                   Folder(".oO[analysis]Oo._September7combination", "Observed", 2, analysis, subdir, plotname="limit_lumi77.45_7813_101,-1.0,1.0_101,-0.02,0.02.root", graphnumber=0, repmap=repmap, linestyle=1, linewidth=2),
                   Folder(".oO[analysis]Oo._September7combination", "Expected", 2, analysis, subdir, plotname="limit_lumi77.45_7813_101,-1.0,1.0_101,-0.02,0.02.root", graphnumber=-1, repmap=repmap, linestyle=7, linewidth=2),
                   Folder(".oO[analysis]Oo._September7combination", "#splitline{Observed,}{2016+2017}", 1, analysis, subdir, plotname="limit_lumi77.45_101,-1.0,1.0_101,-0.02,0.02.root", graphnumber=0, repmap=repmap, linestyle=1, linewidth=1),
                   Folder(".oO[analysis]Oo._September7combination", "#splitline{Expected,}{2016+2017}", 1, analysis, subdir, plotname="limit_lumi77.45_101,-1.0,1.0_101,-0.02,0.02.root", graphnumber=-1, repmap=repmap, linestyle=7, linewidth=1),
                  ]

        mg = ROOT.TMultiGraph("limit", "")
        for folder in folders:
            mg.Add(folder.graph)

        drawlineskwargs = commondrawlineskwargs.copy()
        drawlineskwargs["xpostext"] = CLtextposition
        drawlineskwargs["arbitraryparameter"] = analysis, markerposition

        CLpadindex = 0
        for _ in xdivides:
          if CLtextposition > _:
            CLpadindex += 1
        drawlineskwargs["textsize"] = (
          0.09 / (3*(xboundaries[CLpadindex+1] - xboundaries[CLpadindex]))
        )

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
        logpads[-1].SetRightMargin(rightmargin * 1 / (xboundaries[-1] - xboundaries[-2]))
        logpads[0].SetLeftMargin(leftmargin * 1 / (xboundaries[1] - xboundaries[0]))
        style.subfig(letter, textsize=.11, x1=.79, x2=.83, y1=.80, y2=.84)

        for i, (linearpad, mg) in enumerate(izip(linearpads, mgs)):
            linearpad.cd()
            mg.Draw("al")
            linearpad.SetBottomMargin(bottommargin*2)

        linearpads[-1].SetRightMargin(rightmargin * 1 / (xboundaries[-1] - xboundaries[-2]))
        linearpads[0].SetLeftMargin(leftmargin * 1 / (xboundaries[1] - xboundaries[0]))

        applystyle(mgs, mglogs, folders, xboundaries, xdivides, ydivide)

        for i, (linearpad, mg) in enumerate(izip(linearpads, mgs)):
            linearpad.cd()
            paddrawlineskwargs = drawlineskwargs.copy()
            paddrawlineskwargs["arbitraryparameter"] = paddrawlineskwargs["arbitraryparameter"] + (i,)
            if i != CLpadindex: paddrawlineskwargs["yshift68"] = paddrawlineskwargs["yshift95"] = 100
            drawlines(**paddrawlineskwargs)

        legendpad.cd()
        l = ROOT.TLegend(*legendposition)
        l.SetBorderSize(0)
        l.SetFillStyle(0)
        for folder in folders:
            folder.addtolegend(l)
        if any(folder.secondcolumn is not None for folder in folders):
            assert all(folder.secondcolumn is not None for folder in folders)
            legend.SetNColumns(2)
        l.SetTextSize(.1)
        l.Draw()
        c.cd()
        style.applycanvasstyle(c)
        style.CMS("", lumi=None, lumitext="5.1 fb^{{-1}} (7 TeV) + 19.7 fb^{{-1}} (8 TeV) + {:.1f} fb^{{-1}} (13 TeV)"
                                                .format(sum(_.dataluminosity for _ in config.productionsforcombine)+config.lumi2015),
                      x1=0.007, x2=1.01, #???
                      drawCMS=False, extratextsize=.039)
        style.CMS("", x1=0.09, x2=1.025, y1=.86, y2=.94, CMStextsize=.06)
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