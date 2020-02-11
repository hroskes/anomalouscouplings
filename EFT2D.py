#!/usr/bin/env python

if __name__ == "__main__":
  import argparse
  p = argparse.ArgumentParser()
  p.add_argument("--scaletoCMS", action="store_true")
  args = p.parse_args()

import itertools
import numpy as np
import os
import pprint
import subprocess

import ROOT

from helperstuff import eft, config, stylefunctions as style
from helperstuff.extendedcounter import ExtendedCounter
from helperstuff.utilities import cd, KeepWhileOpenFile, mkdir_p, PlotCopier, TFile

combinecmd = [
  "combine",
  "--robustFit=0",
  "--algo=grid",
  "--floatOtherPOIs=1",
  "--alignEdges=1",
  "--saveSpecifiedFunc=g1,g2,g4,g1prime2,ghg2,ghg4,kappa,kappa_tilde",
  "--saveSpecifiedNuis=everything_but_binbybin",
  "--saveInactivePOI=1",
  "-P", "{coupling1}",
  "-P", "{coupling2}",
  "--X-rtd", "OPTIMIZE_BOUNDS=0",
  "--X-rtd", "TMCSO_AdaptivePseudoAsimov=0",
  "--cminDefaultMinimizerStrategy=0",
  "--X-rtd", "MINIMIZER_MaxCalls=999999999",
  "-t", "-1",
  "--setParameters=g1=1.0,g4=0.0,g2=0.0,g1prime2=0.0",
  "-V",
  "-v", "3",
  "--saveNLL",
  "--setParametersForGrid=g1=1.0",
  "-M", "MultiDimFit",
  "-m", "125",
  "-d", "workspace_lumi3000.00_scang1.root",
  "--setParameterRanges={parameterranges}",
  "--points={npoints}",
  "-n", "{appendname}",
]

def EFT2D(**kwargs):
  couplings = sorted(set(kwargs.pop("couplings")), key="g1 g2 g1prime2 g4".index)
  npoints = kwargs.pop("npoints")
  plotcopier = kwargs.pop("plotcopier")
  scaletoCMS = kwargs.pop("scaletoCMS")
  assert not kwargs, kwargs

  scale = 137.1/3000 if scaletoCMS else 1
  zmax = 50*scale

  repmap = {
    "coupling1": couplings[0],
    "coupling2": couplings[1],
    "parameterranges": ":".join(coupling+"="+{
      "g1": "0.85,1.1",
      "g2": "-0.05,0.05",
      "g4": "-0.1,0.1",
      "g1prime2": "-0.02,0.02",
    }[coupling] for coupling in couplings),
    "appendname": "_2D_{}points_{}_{}".format(npoints, *couplings),
    "npoints": npoints**2
  }
  cmd = [_.format(**repmap) for _ in combinecmd]

  outname = "higgsCombine{appendname}.MultiDimFit.mH125.root".format(**repmap)

  with cd(os.path.join(config.repositorydir, "scans", "cards_fa3fa2fL1_EFT_writeup_width")):
    with KeepWhileOpenFile(outname+".tmp") as kwof:
      if not kwof: return
      if not os.path.exists(outname):
        subprocess.check_call(cmd)

    with TFile(outname) as f:
      t = f.limit
      NLL = ExtendedCounter()

      couplingfunctions = [{
        "g1": lambda t: eft.deltacz(ghz1=t.g1),
        "g2": lambda t: eft.czz(ghz2=t.g2),
        "g4": lambda t: eft.czztilde(ghz4=t.g4),
        "g1prime2": lambda t: eft.czbox(ghz1prime2=t.g1prime2 * 1e4),
      }[_] for _ in couplings]

      for entry in t:
        NLL[tuple(f(entry) for f in couplingfunctions)] = t.deltaNLL
    NLL.zero()

    palettered = np.array([25./256.,246./256.,1,0])
    palettegreen = np.array([121./256.,198./256.,1,153./256.])
    paletteblue = np.array([218./256.,108./256.,1,150./256.])
    paletteposition = np.array([0.,2.3/zmax, 5.99/zmax, 1.0])
    assert len(palettered) == len(palettegreen) == len(paletteblue) == len(paletteposition)
    ncont = 255
    print ROOT.TColor.CreateGradientColorTable(len(palettered), paletteposition, palettered, palettegreen, paletteblue, ncont)
    ROOT.gStyle.SetNumberContours(ncont)

    g = NLL.TGraph2D()

    c = plotcopier.TCanvas("c", "", 8, 30, 800, 800)
    g.Draw("COLZ")

    xtitle, ytitle = ({
      "g1": "#deltac_{z}",
      "g2": "c_{zz}",
      "g4": "#tilde{c}_{zz}",
      "g1prime2": "c_{z#Box}",
    }[_] for _ in couplings)

    xrange, yrange = ({
      "g1": (0.85,1.1),
      "g2": (-0.05,0.05),
      "g4": (-0.1,0.1),
      "g1prime2": (-0.02,0.02),
    }[coupling] for coupling in couplings)

    h = g.GetHistogram()
    h.Scale(scale)

    h.GetXaxis().SetTitle(xtitle)
    #g.GetXaxis().SetRangeUser(*xrange)
    h.GetYaxis().SetTitle(ytitle)
    #g.GetYaxis().SetRangeUser(*yrange)
    h.GetZaxis().SetTitle("#minus2 #Deltaln L")
    h.GetZaxis().SetTitleOffset(1.4)
    h.GetZaxis().SetRangeUser(0, zmax)

    h2 = h.Clone()
    h2.SetContour(2)
    h2.SetContourLevel(0, 2.295748928898637)
    h2.SetContourLevel(1, 5.991464547107987)
    h2.SetLineColor(1)
    h2.Draw("cont2 SAME")

    style.applycanvasstyle(c)
    c.SetRightMargin(0.15)
    style.applyaxesstyle(g)

    SM = ROOT.TGraph(1, np.array([0 if couplings[0] != "g1" else 1]), np.array([0]))
    SM.SetMarkerStyle(30)
    SM.SetMarkerSize(2)
    SM.Draw("P")

    mkdir_p(os.path.join(config.plotsbasedir, "limits", "fa3fa2fL1_EFT_writeup_width", "scaletoCMS"))

    plotname = os.path.join(config.plotsbasedir, "limits", "fa3fa2fL1_EFT_writeup_width", "scaletoCMS" if args.scaletoCMS else "", "2D_{coupling1}_{coupling2}{append}.{ext}")

    for repmap["ext"] in "png root pdf C".split():
      c.SaveAs(plotname.format(append="_nolegend", **repmap))

    dummy1 = ROOT.TGraph(1, np.array([1]), np.array([1]))
    dummy1.SetLineStyle(5)
    dummy2 = ROOT.TGraph(1, np.array([1]), np.array([1]))
    dummy2.SetLineStyle(7)

    l = ROOT.TLegend(.35, .75, .65, .9)
    l.SetBorderSize(0)
    l.SetFillStyle(0)
    l.AddEntry(SM, "SM", "p")
    l.AddEntry(dummy1, "68% CL", "l")
    l.AddEntry(dummy2, "95% CL", "l")
    l.Draw()

    for repmap["ext"] in "png root pdf C".split():
      c.SaveAs(plotname.format(append="", **repmap))

    if list(couplings) != ["g1", "g4"]: return

    c = plotcopier.TCanvas()
    l.Draw()
    for ext in "png pdf root C".split():
      c.SaveAs(os.path.join(config.plotsbasedir, "limits", "fa3fa2fL1_EFT_writeup_width", "scaletoCMS" if args.scaletoCMS else "", "2D_legend."+ext))

if __name__ == "__main__":
  with PlotCopier() as pc:
    for couplings in itertools.combinations("g1 g2 g4 g1prime2".split(), 2):
      EFT2D(plotcopier=pc, npoints=51, couplings=couplings, **args.__dict__)
