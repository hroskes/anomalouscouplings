#!/usr/bin/env python

import argparse

if __name__ == "__main__":
  p = argparse.ArgumentParser()
  p.add_argument("fai", choices="fa3 fa2 fL1 fL1Zg".split())
  args = p.parse_args()

import contextlib, itertools, os
import numpy as np

import ROOT

from TemplateBuilder.TemplateBuilder.fileio import RootFile, RootFiles

import helperstuff.stylefunctions as style

from helperstuff.config import plotsbasedir
from helperstuff.utilities import PlotCopier, reglob

def allthesame(iterable):
  s = set(iterable)
  assert len(s) == 1
  return s.pop()

@contextlib.contextmanager
def RootFilesOrDummies(*filenames, **kwargs):
  if not filenames: yield []; return

  commonargs = kwargs.pop("commonargs", ())
  if kwargs: raise ValueError("Unknown kwargs: "+" ".join(kwargs))

  if filenames[0] is None:
    with RootFilesOrDummies(*filenames[1:], commonargs=commonargs) as morefs:
      yield [None] + morefs
  else:
    with RootFile(filenames[0], *commonargs) as f, RootFilesOrDummies(*filenames[1:], commonargs=commonargs) as morefs:
      yield [f]+morefs

def mergeidenticalscans(outfile, *infiles):
  othercouplingfiles = {
    coupling: [
      infile.replace("limit_", coupling+"_") if "fixothers" not in infile else None
      for infile in infiles
      if all(
        os.path.exists(infile.replace("limit_", coupling+"_"))
        for infile in infiles
        if "fixothers" not in infile
      )
    ] for coupling in ("fa3", "fa2", "fL1", "fL1Zg", "fa1")
  }

  othercouplingfs = {}
  with \
    RootFiles(*infiles) as fs, \
    RootFilesOrDummies(*othercouplingfiles["fa1"]) as othercouplingfs["fa1"], \
    RootFilesOrDummies(*othercouplingfiles["fa2"]) as othercouplingfs["fa2"], \
    RootFilesOrDummies(*othercouplingfiles["fa3"]) as othercouplingfs["fa3"], \
    RootFilesOrDummies(*othercouplingfiles["fL1"]) as othercouplingfs["fL1"], \
    RootFilesOrDummies(*othercouplingfiles["fL1Zg"]) as othercouplingfs["fL1Zg"] \
  :
    for k, v in othercouplingfs.copy().iteritems():
      if not v: del othercouplingfs[k]

    cs = [f.c1 for f in fs]
    frames, multigraphs, legends, text1, text2, line1, line2, text3, text4 = itertools.izip_longest(*(c.GetListOfPrimitives() for c in cs))
    listsofgraphs = [mg.GetListOfGraphs() for mg in multigraphs]
    assert all(len(_) == 1 for _ in listsofgraphs)
    graphs = [_[0] for _ in listsofgraphs]

    _, xxs, yys = itertools.izip(*(itertools.izip(*itertools.izip(xrange(g.GetN()), g.GetX(), g.GetY())) for g in graphs))

    othercouplingcs = {k: [f.c1 if f else None for f in v] for k, v in othercouplingfs.iteritems()}
    othercouplingmultigraphs = {k: [c.GetListOfPrimitives()[1] if c else None for c in v] for k, v in othercouplingcs.iteritems()}
    othercouplinglistsofgraphs = {k: [mg.GetListOfGraphs() if mg else None for mg in v] for k, v in othercouplingmultigraphs.iteritems()}
    assert all(len(_) == 1 for v in othercouplinglistsofgraphs.itervalues() for _ in v if _)
    othercouplinggraphs = {k: [_[0] if _ else None for _ in v] for k, v in othercouplinglistsofgraphs.iteritems()}
    othercouplingxxs = {k: zip(*(itertools.izip(*itertools.izip(xrange(g.GetN()), g.GetX(), g.GetY())) if g else [None]*3 for g in v))[1] for k, v in othercouplinggraphs.iteritems()}
    othercouplingyys = {k: zip(*(itertools.izip(*itertools.izip(xrange(g.GetN()), g.GetX(), g.GetY())) if g else [None]*3 for g in v))[2] for k, v in othercouplinggraphs.iteritems()}

    for fai, faixxs in othercouplingxxs.iteritems():
      for limitxxs, faixxs in itertools.izip(xxs, faixxs):
        if limitxxs != faixxs is not None:
          raise ValueError("Not the same:\n{}\n{}".format(limitxxs, faixxs))
    yyswithfais = [
      [
        (yys[i][j], {k: othercouplingyys[k][i][j] if othercouplingyys[k][i] else 0 for k in othercouplingyys})
        for j in range(len(yys[i]))
      ] for i in range(len(yys))
    ]

    newxs = np.array(sorted(set.union(*(set(_) for _ in xxs))))
    newyswithfais = [min(y for xx, yy in itertools.izip(xxs, yyswithfais) for x, y in itertools.izip(xx, yy) if x == target) for target in newxs]
    newys = np.array([y[0] for y in newyswithfais])
    newfais = {k: np.array([y[1][k] for y in newyswithfais]) for k in othercouplingfs}
    newn = len(newxs)

    for x, y in itertools.izip_longest(newxs, newys): print x, y

    newg = ROOT.TGraph(newn, newxs, newys)
    faigs = {k: ROOT.TGraph(newn, newxs, faiys) for k, faiys in newfais.iteritems()}

    for _ in [newg] + faigs.values():
      _.SetLineStyle(allthesame(g.GetLineStyle() for g in graphs))
      _.SetLineColor(allthesame(g.GetLineColor() for g in graphs))
      _.SetLineWidth(allthesame(g.GetLineWidth() for g in graphs))
      _.SetMarkerStyle(allthesame(g.GetMarkerStyle() for g in graphs))
      _.SetMarkerColor(allthesame(g.GetMarkerColor() for g in graphs))

    newmg = ROOT.TMultiGraph()
    faimgs = {k: ROOT.TMultiGraph() for k in faigs}
    newmg.Add(newg)
    for k in faigs:
      faimgs[k].Add(faigs[k])

    c = pc.TCanvas("c1", "", 8, 30, 800, 800)
    style.applycanvasstyle(c)

    newmg.Draw("AL")
    style.applyaxesstyle(newmg)
    newmg.GetXaxis().SetRangeUser(-1, 1)

    todraw = []
    for _ in legends, text1, text2, line1, line2, text3, text4:
      _ = [thing for thing in _ if thing][0]
      _.Draw()
      _.SetBit(ROOT.kCanDelete, False)

    c.SaveAs(outfile+".png")
    c.SaveAs(outfile+".root")
    c.SaveAs(outfile+".pdf")
    c.SaveAs(outfile+".C")

    for k, faimg in faimgs.iteritems():
      faimg.Draw("AL")
      style.applyaxesstyle(faimg)
      faimg.GetXaxis().SetRangeUser(-1, 1)
      for _ in legends, text1, text2:
        [thing for thing in _ if thing][0].Draw()

      c.SaveAs(outfile.replace("limit_", k+"_")+".png")
      c.SaveAs(outfile.replace("limit_", k+"_")+".root")
      c.SaveAs(outfile.replace("limit_", k+"_")+".pdf")
      c.SaveAs(outfile.replace("limit_", k+"_")+".C")

if __name__ == "__main__":
  pc = PlotCopier()

  with pc:
    mergeidenticalscans(
      os.path.join(plotsbasedir, "limits/fa3fa2fL1fL1Zg_decay_fixsign/limit_lumi300.00_Untagged_scan"+args.fai+"_merged"),
      *reglob(
         os.path.join(plotsbasedir, "limits/fa3fa2fL1fL1Zg_decay_fixsign/"),
         "limit_lumi300.00_Untagged_scan"+args.fai+"(_(f(a1|a3|a2|L1|L1Zg),){4}(f(a1|a3|a2|L1|L1Zg))|_fixothers|).root",
      )
    )
