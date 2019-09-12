#!/usr/bin/env python

import argparse

if __name__ == "__main__":
  p = argparse.ArgumentParser()
  p.add_argument("fai", choices="fa3 fa2 fL1 fL1Zg".split())
  g = p.add_mutually_exclusive_group(required=True)
  g.add_argument("--decay", action="store_true")
  g.add_argument("--boosted", action="store_true")
  g.add_argument("--STXS", action="store_true")
  g.add_argument("--morecategories", action="store_true")
  g.add_argument("--shapesystematics", action="store_true")
  g.add_argument("--newyields", action="store_true")
  g.add_argument("--applySIPcut", action="store_true")
  g.add_argument("--applySIPcut-untaggednosysts", action="store_true")
  p.add_argument("--zoom", action="store_true")
  args = p.parse_args()

import contextlib, itertools, os
import numpy as np

import ROOT

from TemplateBuilder.TemplateBuilder.fileio import RootFile, RootFiles

import helperstuff.stylefunctions as style

from helperstuff.config import plotsbasedir
from helperstuff.plotlimits import drawlines
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

def mergeidenticalscans(outfile, *infiles, **kwargs):
  ngraphs = kwargs.pop("ngraphs", 1)
  zoom = kwargs.pop("zoom", False)
  assert not kwargs, kwargs

  nfiles = len(infiles)

  xtitle = None
  ytitle = None
  othercouplingytitle = None

  othercouplings = [coupling for coupling in ("fa3", "fa2", "fL1", "fL1Zg", "fa1") if "scan"+coupling+"_" not in outfile]

  newgs = []
  faigss = []
  indices = []
  badfiles = []

  for graphindex in range(ngraphs):
    xxs = []
    yyswithfais = []

    for i, infile in enumerate(infiles):
      print i+1, "/", nfiles, os.path.basename(infile)
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
        assert len(listofgraphs) == ngraphs, (len(listofgraphs), ngraphs)

        g = listofgraphs[graphindex]
        _, xx, yy = itertools.izip(*itertools.izip(xrange(g.GetN()), g.GetX(), g.GetY()))
  
        othercouplingc = {k: v.c1 if v else None for k, v in othercouplingf.iteritems()}
        othercouplingmg = {k: v.GetListOfPrimitives()[1] if v else None for k, v in othercouplingc.iteritems()}
        othercouplinglistofgraphs = {k: v.GetListOfGraphs() if v else None for k, v in othercouplingmg.iteritems()}
        assert all(len(v) == ngraphs for v in othercouplinglistofgraphs.itervalues() if v is not None)
        othercouplingg = {k: v[graphindex] if v else None for k, v in othercouplinglistofgraphs.iteritems()}
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
  
        if "fixothers" in infile:
          for x, (y, ywithfais) in itertools.izip_longest(xx, yyswithfais[-1]):
            assert not any(ywithfais.itervalues())
            ywithfais["fa1"] = 1 - abs(x)
  
        if xtitle is None and any(othercouplingmg.itervalues()):
          xtitle = mg.GetXaxis().GetTitle()
          ytitle = mg.GetYaxis().GetTitle()
          othercouplingytitle = {k: v.GetYaxis().GetTitle() for k, v in othercouplingmg.iteritems() if v is not None}
  
    newxs = np.array(sorted(set.union(*(set(_) for _ in xxs))))
    newyswithfais = [min(y for xx, yy in itertools.izip(xxs, yyswithfais) for x, y in itertools.izip(xx, yy) if x == target) for target in newxs]
    newys = np.array([y[0] for y in newyswithfais])
    newfais = {k: np.array([y[1][k] for y in newyswithfais]) for k in othercouplings}
    if "fa3" in newfais: newfais["fa3"] = abs(newfais["fa3"]) #because sign doesn't matter
    newn = len(newxs)
  
    indices.append([{i for i, (xx, yy) in enumerate(itertools.izip(xxs, yyswithfais)) for x, (y, fais) in itertools.izip(xx, yy) if x == target and np.isclose(y, miny)} for target, miny in itertools.izip_longest(newxs, newys)])
  
    fmt = " ".join(["{:>10}"] * (2+len(newfais))) + " {}"
    print fmt.format("x", "y", *[k for k, v in sorted(newfais.iteritems())] + ["file indices"])
    fmt = " ".join(["{:10.3g}"] * (2+len(newfais))) + " {}"
    for xyfaisindices in itertools.izip_longest(newxs, newys, *[v for k, v in sorted(newfais.iteritems())] + [(", ".join(str(idx) for idx in _) for _ in indices[-1])]): print fmt.format(*xyfaisindices)
  
    if 0 not in newxs:
      assert args.newyields
      newys = np.concatenate((newys[newxs < 0], [0], newys[newxs > 0]))
      newfais = {k: np.concatenate((v[newxs < 0], [0], v[newxs > 0])) for k, v in newfais.iteritems()}
      newxs = np.concatenate((newxs[newxs < 0], [0], newxs[newxs > 0]))
      newn += 1
  
    newys -= min(newys)
  
    newg = ROOT.TGraph(newn, newxs, newys)
    faigs = {k: ROOT.TGraph(newn, newxs, faiys) for k, faiys in newfais.iteritems()}

    for _ in [newg] + faigs.values():
      _.SetLineStyle(g.GetLineStyle())
      _.SetLineColor(g.GetLineColor())
      _.SetLineWidth(g.GetLineWidth())
      _.SetMarkerStyle(g.GetMarkerStyle())
      _.SetMarkerColor(g.GetMarkerColor())

    newgs.append(newg)
    faigss.append(faigs)

  newmg = ROOT.TMultiGraph()
  faimgs = {k: ROOT.TMultiGraph() for k in faigs}

  for newg in newgs:
    newmg.Add(newg)
  for faigs in faigss:
    for k in faigs:
      faimgs[k].Add(faigs[k])

  c = pc.TCanvas("c1", "", 8, 30, 800, 800)
  style.applycanvasstyle(c)

  newmg.Draw("AL")
  style.applyaxesstyle(newmg)
  if zoom:
    newmg.GetXaxis().SetRangeUser(-0.02, 0.02)
  else:
    newmg.GetXaxis().SetRangeUser(-1, 1)
  newmg.GetXaxis().SetTitle(xtitle)
  newmg.GetYaxis().SetTitle(ytitle)
  newmg.SetMinimum(0)

  style.CMS("Preliminary", lumi=None, lumitext="{:.1f} fb^{{-1}} (13 TeV)".format(137.1))
  drawlines()

  if zoom: outfile = outfile.replace("_101,-1.0,1.0", "")
  c.SaveAs(outfile+".png")
  c.SaveAs(outfile+".root")
  c.SaveAs(outfile+".pdf")
  c.SaveAs(outfile+".C")

  for k, faimg in faimgs.iteritems():
    faimg.Draw("AL")
    style.applyaxesstyle(faimg)
    faimg.GetXaxis().SetRangeUser(-1, 1)
    faimg.GetXaxis().SetTitle(xtitle)
    print othercouplingytitle
    faimg.GetYaxis().SetTitle(othercouplingytitle[k])

    style.CMS("Preliminary", lumi=None, lumitext="{:.1f} fb^{{-1}} (13 TeV)".format(137.1))

    c.SaveAs(outfile.replace("limit_", k+"_")+".png")
    c.SaveAs(outfile.replace("limit_", k+"_")+".root")
    c.SaveAs(outfile.replace("limit_", k+"_")+".pdf")
    c.SaveAs(outfile.replace("limit_", k+"_")+".C")

  allindices = sum(indices, [])
  allindices.sort(key=lambda x: len(x))
  neededindices = set()
  for indicesatpoint in allindices:
    if not indicesatpoint.intersection(neededindices):
      neededindices.add(indicesatpoint.pop())
  print
  print "The following scans are needed for the final result:"
  print neededindices
  neededfiles = sorted(infiles[idx] for idx in neededindices)
  for _ in neededfiles:
    print "  ", _
  print

if __name__ == "__main__":
  pc = PlotCopier()
  functionkwargs = {"zoom": args.zoom}

  if args.decay:
    functionargs = [
      os.path.join(plotsbasedir, "limits/fa3fa2fL1fL1Zg_decay_fixsign/limit_lumi300.00_Untagged_scan"+args.fai+"_merged"),
    ] + sorted(
      reglob(
         os.path.join(plotsbasedir, "limits/fa3fa2fL1fL1Zg_decay_fixsign/"),
         "limit_lumi300.00_Untagged_scan"+args.fai+"(_(f(a1|a3|a2|L1|L1Zg),){4}(f(a1|a3|a2|L1|L1Zg))|_fixothers|)(_(CMS_zz4l_fai?[0-9]_relative=[0-9.-]*,?)*)?.root",
      ) + reglob(
         os.path.join(plotsbasedir, "limits/fa3fa2fL1fL1Zg_decay_fixsign/gridscan/"),
         "limit_lumi300.00_Untagged_scan"+args.fai+"(_(f(a1|a3|a2|L1|L1Zg),){4}(f(a1|a3|a2|L1|L1Zg))|_fixothers|)(_(CMS_zz4l_fai?[0-9]_relative=[0-9.-]*,?)*)?.root",
      )
    )
  else:
    scanrangeformatch = None
    plotnameformatch = None
    insertinmiddle = ""
    if args.decay:
      folder = "fa3fa2fL1fL1Zg_yieldsystematics"
      plotname = "limit_lumi137.10_scan"+args.fai
      scanrange = ""
    if args.boosted:
      folder = "fa3fa2fL1fL1Zg_boosted_yieldsystematics"
      plotname = "limit_lumi137.10_scan"+args.fai
      scanrange = "_101,-1.0,1.0_101,-0.02,0.02"
    if args.STXS:
      folder = "fa3fa2fL1fL1Zg_STXS_yieldsystematics"
      plotname = "limit_lumi137.10_scan"+args.fai
      scanrange = "_101,-1.0,1.0_101,-0.02,0.02"
    if args.morecategories:
      folder = "fa3fa2fL1fL1Zg_morecategories_yieldsystematics"
      plotname = "limit_lumi137.10_scan"+args.fai
      scanrange = "_101,-1.0,1.0_101,-0.02,0.02"
    if args.shapesystematics:
      folder = "fa3fa2fL1fL1Zg_morecategories_shapesystematics"
      plotname = "limit_lumi137.10_scan"+args.fai
      scanrange = "_101,-1.0,1.0_101,-0.02,0.02"
    if args.newyields:
      folder = "fa3fa2fL1fL1Zg_morecategories_newyields"
      plotname = "limit_lumi137.10_scan"+args.fai
      scanrange = "_26,-1.0,1.0_26,-0.02,0.02"
    if args.applySIPcut:
      folder = "fa3fa2fL1fL1Zg_morecategories_applySIPcut"
      plotname = "limit_lumi137.10_scan"+args.fai
      scanrange = "_101,-1.0,1.0_101,-0.02,0.02"
      scanrangeformatch = "(_101,-0.02,0.02)?"
      functionkwargs["ngraphs"] = 2
    if args.applySIPcut_untaggednosysts:
      folder = "fa3fa2fL1fL1Zg_morecategories_applySIPcut"
      plotname = "limit_lumi137.10_Untagged_scan"+args.fai+"_nosystematics"
      plotnameformatch = "limit_lumi137.10_Untagged_scan"+args.fai+"(?=.*_nosystematics)(_nosystematics)?"
      scanrange = "_101,-1.0,1.0_101,-0.02,0.02"
      scanrangeformatch = "(_nosystematics)?(_101,-0.02,0.02)?"
      insertinmiddle = "(_nosystematics)?"
      functionkwargs["ngraphs"] = 2
    if scanrangeformatch is None: scanrangeformatch = scanrange
    if plotnameformatch is None: plotnameformatch = plotname
    if insertinmiddle is None: insertinmiddle = ""

    functionargs = [
      os.path.join(plotsbasedir, "limits", folder, plotname+scanrange+"_merged"),
    ] + sorted(
      reglob(
         os.path.join(plotsbasedir, "limits", folder, ""),
         plotnameformatch+"(_(f(a1|a3|a2|L1|L1Zg),){4}(f(a1|a3|a2|L1|L1Zg))|_fixothers|)"+insertinmiddle+"(_(CMS_zz4l_fai?[0-9]_relative=[0-9.-]*,?)*)?"+scanrangeformatch+".root",
         hastomatch=True,
      ) + reglob(
         os.path.join(plotsbasedir, "limits", folder, "gridscan/"),
         plotnameformatch+"(_(f(a1|a3|a2|L1|L1Zg),){4}(f(a1|a3|a2|L1|L1Zg))|_fixothers|)"+insertinmiddle+"(_(CMS_zz4l_fai?[0-9]_relative=[0-9.-]*,?)*)?"+scanrangeformatch+".root",
         okifnofolder=True,
      )
    )

  with pc:
    mergeidenticalscans(*functionargs, **functionkwargs)
