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
  nfiles = len(infiles)

  xxs = []
  yyswithfais = []

  othercouplings = [coupling for coupling in ("fa3", "fa2", "fL1", "fL1Zg", "fa1") if "scan"+coupling not in outfile]

  fileforotherstuff = None

  for i, infile in enumerate(infiles, start=1):
    if i % 100 == 0 or i == nfiles: print i, "/", nfiles
    othercouplingfile = {
      coupling: infile.replace("limit_", coupling+"_") if "fixothers" not in infile else None
      for coupling in othercouplings
    }
    othercouplingf = {}
    with RootFilesOrDummies(infile, *(othercouplingfile.get(coupling, None) for coupling in ("fa3", "fa2", "fL1", "fL1Zg", "fa1"))) as (
      f, othercouplingf["fa3"], othercouplingf["fa2"], othercouplingf["fL1"], othercouplingf["fL1Zg"], othercouplingf["fa1"],
    ):
      c = f.c1
      mg = c.GetListOfPrimitives()[1]

      listofgraphs = mg.GetListOfGraphs()
      assert len(listofgraphs) == 1
      g = listofgraphs[0]
      _, xx, yy = itertools.izip(*itertools.izip(xrange(g.GetN()), g.GetX(), g.GetY()))

      othercouplingc = {k: v.c1 if v else None for k, v in othercouplingf.iteritems()}
      othercouplingmg = {k: v.GetListOfPrimitives()[1] if v else None for k, v in othercouplingc.iteritems()}
      othercouplinglistofgraphs = {k: v.GetListOfGraphs() if v else None for k, v in othercouplingmg.iteritems()}
      assert all(len(v) == 1 for v in othercouplinglistofgraphs.itervalues() if v is not None)
      othercouplingg = {k: v[0] if v else None for k, v in othercouplinglistofgraphs.iteritems()}
      othercouplingxx = {k: zip(*itertools.izip(xrange(v.GetN()), v.GetX()))[1] if v else None for k, v in othercouplingg.iteritems()}
      othercouplingyy = {k: zip(*itertools.izip(xrange(v.GetN()), v.GetY()))[1] if v else None for k, v in othercouplingg.iteritems()}

      for fai, faixx in othercouplingxx.iteritems():
        if xx != faixx is not None:
          raise ValueError("Not the same:\n{}\n{}".format(xx, faixx))

      xxs.append(xx)
      yyswithfais.append([
        (yy[j], {k: othercouplingyy[k][j] if othercouplingyy[k] else 0 for k in othercouplingyy})
        for j in range(len(yy))
      ])

      if fileforotherstuff is None:
        if len(c.GetListOfPrimitives()) == 9:
          fileforotherstuff = infile

  newxs = np.array(sorted(set.union(*(set(_) for _ in xxs))))
  newyswithfais = [min(y for xx, yy in itertools.izip(xxs, yyswithfais) for x, y in itertools.izip(xx, yy) if x == target) for target in newxs]
  newys = np.array([y[0] for y in newyswithfais])
  newfais = {k: np.array([y[1][k] for y in newyswithfais]) for k in othercouplings}
  newn = len(newxs)

  for x, y in itertools.izip_longest(newxs, newys): print x, y

  newg = ROOT.TGraph(newn, newxs, newys)
  faigs = {k: ROOT.TGraph(newn, newxs, faiys) for k, faiys in newfais.iteritems()}

  for _ in [newg] + faigs.values():
    _.SetLineStyle(g.GetLineStyle())
    _.SetLineColor(g.GetLineColor())
    _.SetLineWidth(g.GetLineWidth())
    _.SetMarkerStyle(g.GetMarkerStyle())
    _.SetMarkerColor(g.GetMarkerColor())

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
  newmg.SetMinimum(0)

  with RootFilesOrDummies(fileforotherstuff) as (f,):
    if f is not None:
      _, _, legend, text1, text2, line1, line2, text3, text4 = f.c1.GetListOfPrimitives()
      for _ in legend, text1, text2, line1, line2, text3, text4:
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
      if f is not None:
        for _ in legend, text1, text2:
          _.Draw()

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
         "limit_lumi300.00_Untagged_scan"+args.fai+"(_(f(a1|a3|a2|L1|L1Zg),){4}(f(a1|a3|a2|L1|L1Zg))|_fixothers|)(_(CMS_zz4l_fai?[0-9]_relative=[0-9.-]*,?)*)?.root",
      )
    )
