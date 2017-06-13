#!/usr/bin/env python

from array import array
import glob
import os
import pipes
import re
import shutil
import subprocess

import ROOT

from helperstuff import config
from helperstuff.enums import analyses, Analysis
from helperstuff.utilities import mkdir_p, mkdtemp

from projections import Projections

def doanimationstep(plotfile, markerposition, saveas):
  f = ROOT.TFile(plotfile)
  try:
    c = f.c1
  except:
    f.ls()
    raise
  c.cd()
  markerposition = [array('d', [_]) for _ in markerposition]
  marker = ROOT.TGraph(1, *markerposition)
  marker.SetMarkerStyle(20)
  marker.SetMarkerColor(4)
  marker.SetMarkerSize(3)
  marker.Draw("P")
  c.SaveAs(saveas)

class WrapTree(ROOT.TChain):
  @property
  def muV_scaled(self):
    return self.r_VVH
  @property
  def muf_scaled(self):
    return self.r_ffH

class Projections_WIN(Projections):
  @classmethod
  def scantreeforanimations(cls, analysis):
    analysis = Analysis(analysis)
    scantree = WrapTree("limit")
    for filename in glob.glob(os.path.join(config.repositorydir, "CMSSW_7_6_5", "src", "HiggsAnalysis", "HZZ4l_Combination",
                                       "CreateDatacards", "cards_{}_Feb28_mu".format(analysis), "higgsCombine_obs_*.root")):
        if re.match("higgsCombine_obs_lumi[0-9.]*(_[0-9]*,[-.0-9]*,[-.0-9]*)*.MultiDimFit.mH125.root", os.path.basename(filename)):
            scantree.Add(filename)
    return scantree


def animatelimits(analysis):
  tmpdir = mkdtemp()
  animation = Projections_WIN.animationstepsforniceplots(analysis)
  convertcommand = ["gm", "convert", "-loop", "0"]
  lastdelay = None
  plotfile = os.path.join(config.plotsbasedir, "limits", ".bkp", "bkp_FebruaryMarch", "{}_Feb28_mu".format(analysis), "limit_lumi35.8671_100,-1.0,1.0_100,-0.02,0.02.root")
  for i, step in enumerate(animation):
    if step.delay != lastdelay:
      convertcommand += ["-delay", str(step.delay)]
    convertcommand += ["-trim", os.path.join(tmpdir, "{}.pdf".format(i))]
    doanimationstep(
                    plotfile=plotfile,
                    saveas=os.path.join(tmpdir, "{}.pdf".format(i)),
                    markerposition=(step.fai_decay, step.deltaNLL),
                   )

  plotfolder = os.path.join(config.plotsbasedir, "limits", "forWIN", "nocombination")
  mkdir_p(plotfolder)
  finalplot = os.path.join(plotfolder, "{}_limit.gif".format(analysis))
  convertcommand.append(finalplot)
  #http://stackoverflow.com/a/38792806/5228524
  #subprocess.check_call(convertcommand)
  os.system(" ".join(pipes.quote(_) for _ in convertcommand))
  shutil.rmtree(tmpdir)

if __name__ == "__main__":
  for analysis in analyses:
    animatelimits(analysis)
