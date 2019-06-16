#!/usr/bin/env python

import argparse

if __name__ == "__main__":
  p = argparse.ArgumentParser()
  p.add_argument("fai", choices="fa3 fa2 fL1 fL1Zg".split())
  args = p.parse_args()

import itertools, os
import numpy as np

import ROOT

from TemplateBuilder.TemplateBuilder.fileio import RootFiles

import helperstuff.stylefunctions as style

from helperstuff.config import plotsbasedir
from helperstuff.utilities import PlotCopier, reglob

def allthesame(iterable):
  s = set(iterable)
  assert len(s) == 1
  return s.pop()

def mergeidenticalscans(outfile, *infiles):
  with RootFiles(*infiles) as fs:
    cs = [f.c1 for f in fs]
    frames, multigraphs, legends, text1, text2, line1, line2, text3, text4 = itertools.izip_longest(*(c.GetListOfPrimitives() for c in cs))
    listsofgraphs = [mg.GetListOfGraphs() for mg in multigraphs]
    assert all(len(_) == 1 for _ in listsofgraphs)
    graphs = [_[0] for _ in listsofgraphs]

    _, xxs, yys = itertools.izip(*(itertools.izip(*itertools.izip(xrange(g.GetN()), g.GetX(), g.GetY())) for g in graphs))

    newxs = np.array(sorted(set.union(*(set(_) for _ in xxs))))
    newys = np.array([min(y for xx, yy in itertools.izip(xxs, yys) for x, y in itertools.izip(xx, yy) if x == target) for target in newxs])
    newn = len(newxs)

    for x, y in itertools.izip_longest(newxs, newys): print x, y

    newg = ROOT.TGraph(newn, newxs, newys)
    newg.SetLineStyle(allthesame(g.GetLineStyle() for g in graphs))
    newg.SetLineColor(allthesame(g.GetLineColor() for g in graphs))
    newg.SetLineWidth(allthesame(g.GetLineWidth() for g in graphs))
    newg.SetMarkerStyle(allthesame(g.GetMarkerStyle() for g in graphs))
    newg.SetMarkerColor(allthesame(g.GetMarkerColor() for g in graphs))

    newmg = ROOT.TMultiGraph()
    newmg.Add(newg)

    c = pc.TCanvas("c1", "", 8, 30, 800, 800)
    newmg.Draw("AL")

    style.applycanvasstyle(c)
    style.applyaxesstyle(newmg)
    newmg.GetXaxis().SetRangeUser(-1, 1)

    for _ in legends, text1, text2, line1, line2, text3, text4:
      _[0].Draw()

    c.SaveAs(outfile+".png")
    c.SaveAs(outfile+".root")
    c.SaveAs(outfile+".pdf")
    c.SaveAs(outfile+".C")

if __name__ == "__main__":
  pc = PlotCopier()

  with pc:
    mergeidenticalscans(
      os.path.join(plotsbasedir, "limits/fa3fa2fL1fL1Zg_decay_editingcombine/limit_lumi300.00_Untagged_scan"+args.fai+"_merged"),
      *reglob(
         os.path.join(plotsbasedir, "limits/fa3fa2fL1fL1Zg_decay_editingcombine/"),
         "limit_lumi300.00_Untagged_scan"+args.fai+"(_(f(a1|a3|a2|L1|L1Zg),){4}(f(a1|a3|a2|L1|L1Zg))|_fixothers|).root",
      )
    )
