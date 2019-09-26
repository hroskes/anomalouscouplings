#!/usr/bin/env python
if __name__ == "__main__":
    import argparse

    def f(x): result = x.split("="); assert len(result)==2, x; return result

    p = argparse.ArgumentParser()
    p.add_argument("kwargs", type=f, nargs="*")
    p.add_argument("--analysis", choices="fa3 fa2 fL1 fL1Zg".split())
    args = p.parse_args()

    kwargs = {k: v for k, v in args.kwargs}
    if "analysis" in kwargs: args.analysis = kwargs["analysis"]

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
from helperstuff.utilities import cache, cd, deprecate, mkdir_p, mkdtemp, PlotCopier, tfiles

from mergeplots import Folder

analyses = "fa3", "fa2", "fL1", "fL1Zg"
setmax = 1
def getplotname(analysis, comparecategories):
    if comparecategories:
        return "limit_lumi77.40_scan{}_101,-1.0,1.0_101,-0.02,0.02_compare_categories_zoom.root".format(analysis)
    return "limit_lumi77.40_scan{}_101,-1.0,1.0_101,-0.02,0.02_compare_zoom.root".format(analysis)

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
        mg.GetXaxis().SetTitleOffset(.6)
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

def PRL_loglinear(**kwargs):
    commondrawlineskwargs = {
                             "logscale": False,  #the lines are in the linear part
                             "xsize": .2,
                             "ysize": .045,
                             "yshift68": .03,
                             "yshift95": .03,
                            }

    markerposition = kwargs.pop("markerposition", None)
    onlyanalysis = kwargs.pop("analysis", None)
    xdivides = sorted(float(_) for _ in kwargs.pop("xdivides", (-.03, .03)))
    assert len(xdivides) == 2, xdivides
    ydivide = float(kwargs.pop("ydivide", 11))
    saveas = kwargs.pop("saveas", None)
    comparecategories = bool(int(kwargs.pop("comparecategories", False)))

    commondrawlineskwargs.update(kwargs)

    saveasdir = os.path.join(config.plotsbasedir, "limits", "fa3fa2fL1fL1Zg_morecategories_finalforthesis")

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
            CLtextposition=.65
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
        biglegend = True
        #assert abs((leftmargin + rightmargin) - (topmargin + bottommargin)) < 1e-5, (leftmargin + rightmargin, topmargin + bottommargin)
        c.SetLeftMargin(0)
        c.SetRightMargin(0)
        c.SetTopMargin(0)
        c.SetBottomMargin(0)
        xboundaries = [0, leftmargin+(1-leftmargin-rightmargin)/3, leftmargin+(1-leftmargin-rightmargin)*2/3, 1]
        yboundaries = [0, bottommargin+(1-bottommargin-topmargin)/2, 1]
        if biglegend:
            if analysis == "fa2":
                legendposition = .4, .57, .8, .87
            elif analysis == "fL1":
                legendposition = .2, .57, .6, .87
            else:
                legendposition = .3, .57, .7, .87
        else:
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

        assert "finalforthesis" in saveasdir
        if analysis == "fa3":
          removepoints1 = []
          removepoints2 = [-0.28, 0.28]
          removepoints3 = [-0.84, -0.82, -0.78, -0.6, -0.56, -0.5, -0.46, -0.42, 0.0008, 0.0012, 0.0016, 0.002, 0.0024, 0.0028, 0.004, 0.0044, 0.0052, 0.006, 0.0068, 0.0072, 0.0076, 0.0092, 0.0096, 0.0108, 0.0116, 0.014, 0.0148, 0.0152, 0.016, 0.0168, 0.0188, 0.0192]
          removepoints4 = []
        if analysis == "fa2":
          removepoints1 = []
          removepoints2 = [0.6, 0.66, 0.76]
          removepoints3 = [-0.14, -0.08, -0.06]
          removepoints4 = [-0.14, -0.12, -0.06, 0.08, 0.12, 0.14, 0.16, 0.2, 0.24, 0.26]
        if analysis == "fL1":
          removepoints1 = []
          removepoints2 = []
          removepoints3 = [-0.014, -0.0124, -0.0092, -0.0084, -0.0052, -0.0048, 0.018]
          removepoints4 = [0.1, 0.12, 0.16, 0.18, 0.22, 0.26, 0.3]
        if analysis == "fL1Zg":
          removepoints1 = []
          removepoints2 = []
          removepoints3 = []
          removepoints4 = [-0.76, -0.66, -0.56]

        if comparecategories:
          folders = [
            Folder("fa3fa2fL1fL1Zg_finalforthesis/", "Float others 3 categories", 2, analysis, subdir, plotname="limit_lumi77.40_scan.oO[analysis]Oo._101,-1.0,1.0_101,-0.02,0.02_merged.root", graphnumber=0, repmap=repmap, linestyle=7, linewidth=2),
            Folder("fa3fa2fL1fL1Zg_finalforthesis/", "Fix others 3 categories", 4, analysis, subdir, plotname="limit_lumi77.40_scan.oO[analysis]Oo._fixothers_101,-1.0,1.0_101,-0.02,0.02.root", graphnumber=0, repmap=repmap, linestyle=7, linewidth=2),
#            Folder("fa3fa2fL1fL1Zg_boosted_finalforthesis/", "Float others w/ boosted", 2, analysis, subdir, plotname="limit_lumi77.40_scan.oO[analysis]Oo._101,-1.0,1.0_101,-0.02,0.02_merged.root", graphnumber=0, repmap=repmap, linestyle=1, linewidth=2),
#            Folder("fa3fa2fL1fL1Zg_boosted_finalforthesis/", "Fix others w/ boosted", 4, analysis, subdir, plotname="limit_lumi77.40_scan.oO[analysis]Oo._fixothers_101,-1.0,1.0_101,-0.02,0.02.root", graphnumber=0, repmap=repmap, linestyle=1, linewidth=2),
#            Folder("fa3fa2fL1fL1Zg_STXS_finalforthesis/", "Float others STXS", 2, analysis, subdir, plotname="limit_lumi77.40_scan.oO[analysis]Oo._101,-1.0,1.0_101,-0.02,0.02_merged.root", graphnumber=0, repmap=repmap, linestyle=8, linewidth=2),
#            Folder("fa3fa2fL1fL1Zg_STXS_finalforthesis/", "Fix others STXS", 4, analysis, subdir, plotname="limit_lumi77.40_scan.oO[analysis]Oo._fixothers_101,-1.0,1.0_101,-0.02,0.02.root", graphnumber=0, repmap=repmap, linestyle=8, linewidth=2),
            Folder("fa3fa2fL1fL1Zg_morecategories_finalforthesis/", "Float others 6 categories", 2, analysis, subdir, plotname="limit_lumi77.40_scan.oO[analysis]Oo._101,-1.0,1.0_101,-0.02,0.02_merged.root", graphnumber=0, repmap=repmap, linestyle=3, linewidth=2),
            Folder("fa3fa2fL1fL1Zg_morecategories_finalforthesis/", "Fix others 6 categories", 4, analysis, subdir, plotname="limit_lumi77.40_scan.oO[analysis]Oo._fixothers_101,-1.0,1.0_101,-0.02,0.02.root", graphnumber=0, repmap=repmap, linestyle=3, linewidth=2),
          ]
        else:
          folders = [
            Folder("fa3fa2fL1fL1Zg_morecategories_finalforthesis/", "Observed, fix others", 4, analysis, subdir, plotname="limit_lumi77.40_scan.oO[analysis]Oo._fixothers_101,-1.0,1.0_101,-0.02,0.02.root", graphnumber=0, repmap=repmap, linestyle=1, linewidth=2, removepoints=removepoints1),
            Folder("fa3fa2fL1fL1Zg_morecategories_finalforthesis/", "Expected, fix others", 4, analysis, subdir, plotname="limit_lumi77.40_scan.oO[analysis]Oo._fixothers_101,-1.0,1.0_101,-0.02,0.02.root", graphnumber=1, repmap=repmap, linestyle=7, linewidth=2, removepoints=removepoints2),
#            Folder("fa3fa2fL1fL1Zg_morecategories_finalforthesis/", "Observed, float others", 2, analysis, subdir, plotname="limit_lumi77.40_scan.oO[analysis]Oo._101,-1.0,1.0_101,-0.02,0.02_merged.root", graphnumber=0, repmap=repmap, linestyle=1, linewidth=2, removepoints=removepoints3),
#            Folder("fa3fa2fL1fL1Zg_morecategories_finalforthesis/", "Expected, float others", 2, analysis, subdir, plotname="limit_lumi77.40_scan.oO[analysis]Oo._101,-1.0,1.0_101,-0.02,0.02_merged.root", graphnumber=1, repmap=repmap, linestyle=7, linewidth=2, removepoints=removepoints4),
            Folder(".oO[analysis]Oo._August17combination/", "Observed, 18-002", 1, analysis, subdir, plotname="limit_lumi77.45_101,-1.0,1.0_101,-0.02,0.02.root", graphnumber=0, repmap=repmap, linestyle=1, linewidth=2),
            Folder(".oO[analysis]Oo._August17combination/", "Expected, 18-002", 1, analysis, subdir, plotname="limit_lumi77.45_101,-1.0,1.0_101,-0.02,0.02.root", graphnumber=1, repmap=repmap, linestyle=7, linewidth=2),
          ]

        mg = ROOT.TMultiGraph("limit", "")
        for folder in folders:
            mg.Add(folder.graph)

        drawlineskwargs = commondrawlineskwargs.copy()
        drawlineskwargs["xpostext"] = CLtextposition

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
        mglog = mg.Clone(mg.GetName()+"_log")
        if markerposition and marker.GetY()[0] <= 0:
            mg.Add(marker, "P")
        mgs = [mg, mg.Clone(mg.GetName()+"_1"), mg.Clone(mg.GetName()+"_2")]
        mglogs = [mglog, mglog.Clone(mglog.GetName()+"_1"), mglog.Clone(mglog.GetName()+"_2")]
        for logpad, mglog in izip(logpads, mglogs):
            logpad.cd()
            logpad.SetLogy()
            mglog.Draw("al")
            logpad.SetTopMargin(topmargin*2)
        logpads[-1].SetRightMargin(rightmargin * 1 / (xboundaries[-1] - xboundaries[-2]))
        logpads[0].SetLeftMargin(leftmargin * 1 / (xboundaries[1] - xboundaries[0]))
        c.cd()
        style.subfig(letter, textsize=.04, x1=.91, x2=.94, y1=.88, y2=.92)

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
            if i != CLpadindex: paddrawlineskwargs["yshift68"] = paddrawlineskwargs["yshift95"] = 100
            drawlines(**paddrawlineskwargs)

        (c if biglegend else legendpad).cd()
        l = ROOT.TLegend(*legendposition)
        l.SetBorderSize(0)
        l.SetFillStyle(0)
        for folder in folders:
            folder.addtolegend(l)
        if any(folder.secondcolumn is not None for folder in folders):
            assert all(folder.secondcolumn is not None for folder in folders)
            legend.SetNColumns(2)
        l.SetTextSize(.04 if biglegend else .1)
        l.Draw()
        c.cd()
        style.applycanvasstyle(c)
        style.CMS("", lumi=None, lumitext="{:.1f} fb^{{-1}} (13 TeV)"
                      .format(sum(_.dataluminosity for _ in config.productionsforcombine if _.year != 2018)),
                      x1=0.007, x2=1.01, #???
                      drawCMS=False, extratextsize=.039)
        style.CMS("", x1=0.09, x2=1.025, y1=.86, y2=.94, CMStextsize=.06, extratextsize=.039)
        yaxislabel(folders[0].ytitle).Draw()

        mkdir_p(os.path.join(saveasdir, "preliminary"))
        mkdir_p(os.path.join(saveasdir, "workinprogress"))
        plotname = getplotname(analysis, comparecategories)
        assert ".root" in plotname

        if saveas is not None:
            c.SaveAs(saveas)
        else:
            for ext in "png eps root pdf".split():
                c.SaveAs(os.path.join(saveasdir, replaceByMap(plotname.replace("root", ext), repmap)))
            with plotcopier.open(os.path.join(saveasdir, replaceByMap(plotname.replace("root", "txt"), repmap)), "w") as f:
                f.write(" ".join(["python"]+[pipes.quote(_) for _ in sys.argv]))
                f.write("\n\n\n\n\n\ngit info:\n\n")
                f.write(subprocess.check_output(["git", "rev-parse", "HEAD"]))
                f.write("\n")
                f.write(subprocess.check_output(["git", "status"]))
                f.write("\n")
                f.write(subprocess.check_output(["git", "diff"]))

        style.CMS("Preliminary", x1=0.12, x2=1.025, y1=.85, y2=.93, drawCMS=False, CMStextsize=.06, extratextsize=.039)

        if saveas is not None:
            assert false
        else:
            for ext in "png eps root pdf".split():
                c.SaveAs(os.path.join(saveasdir, "preliminary", replaceByMap(plotname.replace("root", ext), repmap)))
            with plotcopier.open(os.path.join(saveasdir, "preliminary", replaceByMap(plotname.replace("root", "txt"), repmap)), "w") as f:
                f.write(" ".join(["python"]+[pipes.quote(_) for _ in sys.argv]))
                f.write("\n\n\n\n\n\ngit info:\n\n")
                f.write(subprocess.check_output(["git", "rev-parse", "HEAD"]))
                f.write("\n")
                f.write(subprocess.check_output(["git", "status"]))
                f.write("\n")
                f.write(subprocess.check_output(["git", "diff"]))

        del c.GetListOfPrimitives()[-1]  #delete Preliminary
        style.CMS("Work in progress", x1=0.12, x2=1.025, y1=.85, y2=.93, drawCMS=False, CMStextsize=.06, extratextsize=.039)

        if saveas is not None:
            assert false
        else:
            for ext in "png eps root pdf".split():
                c.SaveAs(os.path.join(saveasdir, "workinprogress", replaceByMap(plotname.replace("root", ext), repmap)))
            with plotcopier.open(os.path.join(saveasdir, "workinprogress", replaceByMap(plotname.replace("root", "txt"), repmap)), "w") as f:
                f.write(" ".join(["python"]+[pipes.quote(_) for _ in sys.argv]))
                f.write("\n\n\n\n\n\ngit info:\n\n")
                f.write(subprocess.check_output(["git", "rev-parse", "HEAD"]))
                f.write("\n")
                f.write(subprocess.check_output(["git", "status"]))
                f.write("\n")
                f.write(subprocess.check_output(["git", "diff"]))

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
    function = PRL_loglinear
    with PlotCopier() as plotcopier:
        for kwargs["analysis"] in "fa3", "fa2", "fL1", "fL1Zg":
            if kwargs["analysis"] != args.analysis is not None: continue
            function(**kwargs)

