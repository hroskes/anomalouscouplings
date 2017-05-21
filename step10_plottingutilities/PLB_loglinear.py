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

from mergeplots import Folder

analyses = "fa3", "fa2", "fL1", "fL1Zg"

def PRL_loglinear(**kwargs):
    commondrawlineskwargs = {
                             "logscale": False,  #the lines are in the linear part
                             "xsize": .3,
                             "ysize": .1,
                             "textsize": .1,
                             "yshift68": .08,
                             "yshift95": -.1,
                            }
    baseplotname = "limit_lumi35.8671.root"
    ydivide = 4.5
    CLtextposition = -0.5
    for kw, kwarg in kwargs.iteritems():
        if kw == "CLtextposition":
            CLtextposition = kwarg
        elif kw == "ydivide":
            ydivide = float(kwarg)
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
        if analysis == "fa3":
            x1legend = .32
        elif analysis == "fa2":
            x1legend = .5
        elif analysis in ("fL1", "fL1Zg"):
            x1legend = .1
        else:
            assert False
        legendposition = x1legend, .3, x1legend+.45, .8

        c = ROOT.TCanvas("c1", "", 8, 30, 1600, 1600)
        leftmargin = .08
        rightmargin = .04 #apply to the individual pads or 1 of the x axis gets cut off
        topmargin = .02
        bottommargin = .1
        assert abs((leftmargin + rightmargin) - (topmargin + bottommargin)) < 1e-5, (leftmargin + rightmargin, topmargin + bottommargin)
        c.SetLeftMargin(leftmargin)
        c.SetRightMargin(0)
        c.Divide(1, 2, 0, 0)
        linearpad = c.cd(2)
        logpad = legendpad = c.cd(1)
        linearpad.SetTicks()
        logpad.SetTicks()

        analysis = Analysis(analysis)
        repmap = {"analysis": str(analysis)}
        subdir = ""
        folders = [
                   Folder(".oO[analysis]Oo._allsysts", "Observed", 2, analysis, subdir, plotname="limit_lumi35.8671_7813_100,-1.0,1.0_100,-0.02,0.02.root", graphnumber=0, repmap=repmap, linestyle=1, linewidth=2),
                   Folder(".oO[analysis]Oo._allsysts", "Expected", 2, analysis, subdir, plotname="limit_lumi35.8671_7813_100,-1.0,1.0_100,-0.02,0.02.root", graphnumber=1, repmap=repmap, linestyle=7, linewidth=2),
                   Folder(".oO[analysis]Oo._allsysts", "Observed, 13 TeV", 1, analysis, subdir, plotname="limit_lumi35.8671_13_100,-1.0,1.0_100,-0.02,0.02.root", graphnumber=0, repmap=repmap, linestyle=1, linewidth=1),
                   Folder(".oO[analysis]Oo._allsysts", "Expected, 13 TeV", 1, analysis, subdir, plotname="limit_lumi35.8671_13_100,-1.0,1.0_100,-0.02,0.02.root", graphnumber=1, repmap=repmap, linestyle=7, linewidth=1),
                  ]

        mg = ROOT.TMultiGraph("limit", "")
        for folder in folders:
            mg.Add(folder.graph)

        setmax = 1

        mglog = mg.Clone()
        logpad.cd()
        logpad.SetLogy()
        mglog.Draw("al")
        mglog.GetXaxis().SetTitle(folders[0].xtitle)
        mglog.GetXaxis().SetRangeUser(-setmax, setmax)
        mglog.GetXaxis().CenterTitle()
        mglog.GetYaxis().CenterTitle()
        mglog.SetMinimum(ydivide)
        mglog.SetMaximum(120)
        style.applyaxesstyle(mglog)
        mglog.GetXaxis().SetLabelSize(.07)
        mglog.GetYaxis().SetLabelSize(.07)
        mglog.GetXaxis().SetTitleSize(.08)
        mglog.GetYaxis().SetTitleSize(.07)
        logpad.SetRightMargin(rightmargin)
        logpad.SetTopMargin(topmargin*2)
        style.subfig(letter, textsize=.08, x1=.89, x2=.93, y1=.89, y2=.93)

        linearpad.cd()
        mg.Draw("al")
        mg.GetXaxis().SetTitle(folders[0].xtitle)
        mg.GetXaxis().SetRangeUser(-setmax, setmax)
        mg.GetXaxis().CenterTitle()
        mg.GetYaxis().CenterTitle()
        mg.SetMinimum(0)
        mg.SetMaximum(ydivide)
        style.applyaxesstyle(mg)
        mg.GetXaxis().SetLabelSize(.07)
        mg.GetYaxis().SetLabelSize(.07)
        mg.GetXaxis().SetTitleSize(.08)
        mg.GetYaxis().SetTitleSize(.07)
        linearpad.SetRightMargin(rightmargin)
        linearpad.SetBottomMargin(bottommargin*2)

        drawlineskwargs = commondrawlineskwargs.copy()
        drawlineskwargs["xpostext"] = CLtextposition
        print drawlineskwargs
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
                                                .format(config.productionforcombine.dataluminosity+config.lumi2015),
                      x1=0.007, x2=1.01, #???
                      drawCMS=False, extratextsize=.03)
        style.CMS("", x1=0.09, x2=1.025, y1=.86, y2=.94, CMStextsize=.05)
        yaxislabel(folders[0].ytitle).Draw()

        saveasdir = os.path.join(config.plotsbasedir, "limits")
        try:
            os.makedirs(saveasdir)
        except OSError:
            pass
        plotname = "{}_{}".format(analysis, baseplotname)
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


@cache
def yaxislabel(label, textsize=.05):
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
    PRL_loglinear(*args, **kwargs)