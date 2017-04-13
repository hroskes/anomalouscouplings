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
    drawlineskwargs = {
                       "PRL": True,
                       "logscale": False,  #the lines are in the linear part
                      }
    legendposition = (.2, .7, .6, .9)
    plotname = "limit_lumi35.8671.root"
    ydivide = 4.5
    for kw, kwarg in kwargs.iteritems():
        if kw == "CLtextposition":
            drawlineskwargs["xpostext"] = kwarg
        elif kw == "legendposition":
            try:
                legendposition = [float(a) for a in kwarg.split(",")]
                if len(legendposition) != 4: raise ValueError
            except ValueError:
                raise ValueError("legendposition has to contain 4 floats separated by commas!")
        elif kw == "ydivide":
            ydivide = float(kwarg)
        else:
            drawlineskwargs[kw] = kwarg

    for k, v in drawlineskwargs.items():
        if k == "xpostext":
            try:
                drawlineskwargs[k] = float(v)
            except ValueError:
                pass
        elif k in ("xmin", "xmax"):
            drawlineskwargs[k] = float(v)

    c = ROOT.TCanvas("c1", "", 8, 30, len(analyses)*1600, 1600)
    c.Divide(len(analyses), 2, 0, 0)
    mgs = {}
    mgslog = {}
    for i, analysis in reversed(list(enumerate(analyses, start=1))):
        linearpad = c.cd(i+len(analyses))
        logpad = c.cd(i)

        analysis = Analysis(analysis)
        repmap = {"analysis": str(analysis)}
        subdir = ""
        folders = [
                   Folder(".oO[analysis]Oo._allsysts", "Observed", 4, analysis, subdir, plotname="limit_lumi35.8671_7813_100,-1.0,1.0_100,-0.02,0.02.root", graphnumber=0, repmap=repmap, linestyle=1, linewidth=2),
                   Folder(".oO[analysis]Oo._allsysts", "Expected", 4, analysis, subdir, plotname="limit_lumi35.8671_7813_100,-1.0,1.0_100,-0.02,0.02.root", graphnumber=1, repmap=repmap, linestyle=7, linewidth=2),
                   Folder(".oO[analysis]Oo._allsysts", "Observed, 13 TeV", 1, analysis, subdir, plotname="limit_lumi35.8671_13_100,-1.0,1.0_100,-0.02,0.02.root", graphnumber=0, repmap=repmap, linestyle=1, linewidth=1),
                   Folder(".oO[analysis]Oo._allsysts", "Expected, 13 TeV", 1, analysis, subdir, plotname="limit_lumi35.8671_13_100,-1.0,1.0_100,-0.02,0.02.root", graphnumber=1, repmap=repmap, linestyle=7, linewidth=1),
                  ]

        mg = mgs[analysis] = ROOT.TMultiGraph("limit", "")
        for folder in folders:
            mg.Add(folder.graph)

        mglog = mgslog[analysis] = mg.Clone()
        logpad.cd()
        logpad.SetLogy()
        mglog.Draw("al")
        mglog.GetXaxis().SetTitle(folders[0].xtitle)
        mglog.GetYaxis().SetTitle(folders[0].ytitle)
        mglog.GetXaxis().SetRangeUser(-1, 1)
        mglog.GetXaxis().CenterTitle()
        mglog.GetYaxis().CenterTitle()
        mglog.SetMinimum(ydivide)
        mglog.SetMaximum(120)
        style.applyaxesstyle(mglog)
        mglog.GetXaxis().SetLabelSize(.08)
        mglog.GetYaxis().SetLabelSize(.08)
        mglog.GetXaxis().SetTitleSize(.08)
        mglog.GetYaxis().SetTitleSize(.08)
        logpad.SetTopMargin(.1)

        linearpad.cd()
        mg.Draw("al")
        mg.GetXaxis().SetTitle(folders[0].xtitle)
        mg.GetYaxis().SetTitle(folders[0].ytitle)
        mg.GetXaxis().SetRangeUser(-1, 1)
        mg.GetXaxis().CenterTitle()
        mg.GetYaxis().CenterTitle()
        mg.SetMinimum(0)
        mg.SetMaximum(ydivide)
        style.applyaxesstyle(mg)
        mg.GetXaxis().SetLabelSize(.08)
        mg.GetYaxis().SetLabelSize(.08)
        mg.GetXaxis().SetTitleSize(.09)
        mg.GetYaxis().SetTitleSize(.09)
        linearpad.SetBottomMargin(.2)

        drawlines(**drawlineskwargs)

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
                  x1=0.03)


    saveasdir = os.path.join(config.plotsbasedir, "limits")
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
    PRL_loglinear(*args, **kwargs)
