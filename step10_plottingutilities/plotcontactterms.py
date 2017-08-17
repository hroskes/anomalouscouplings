#!/usr/bin/env python

from array import array
import os
import tempfile

import ROOT

from helperstuff import config
from helperstuff.copyplots import copyplots
from helperstuff.plotlimits import plotlimits

def plotcontactterms(folder):
  f = ROOT.TFile(os.path.join(config.plotsbasedir, "limits", folder, "limit_lumi300.0_2e2mu_Untagged_10000,-1.0,1.0.root"))
  c1 = f.c1
  g = c1.GetListOfPrimitives()[1]
  left = c1.GetListOfPrimitives()[2].GetListOfGraphs()[0]
  assert left.GetLineColor() == 2
  right = c1.GetListOfPrimitives()[2].GetListOfGraphs()[2]
  assert right.GetLineColor() == 3

  filename = maketree(g, left)
  print plotlimits(os.path.join(config.plotsbasedir, "limits", folder, "limit_feL.root"), "fL1", 0, infilename=filename, drawCMS=False, CMStext="", luminosity=300, xtitle="f_{#epsilonL}cos(#phi_{#lower[-0.25]{#epsilonL}})", scanranges=[(100,-1,1)], killpoints=(-.35, -.3))

  filename = maketree(g, right)
  print plotlimits(os.path.join(config.plotsbasedir, "limits", folder, "limit_feR.root"), "fL1", 0, infilename=filename, drawCMS=False, CMStext="", luminosity=300, xtitle="f_{#epsilonR}cos(#phi_{#lower[-0.25]{#epsilonR}})", scanranges=[(100,-1,1)], killpoints=(.3, .35))

  copyplots(os.path.join("limits", folder))

def maketree(g, line):
  tmpdir = tempfile.mkdtemp()
  filename = os.path.join(tmpdir, "tmp.root")
  f = ROOT.TFile(filename, "CREATE")
  t = ROOT.TTree("limit", "limit")
  CMS_zz4l_fai1 = array("d", [0])
  deltaNLL = array("d", [0])
  nll0 = array("d", [0])
  nll = array("d", [0])
  t.Branch("CMS_zz4l_fai1", CMS_zz4l_fai1, "CMS_zz4l_fai1/D")
  t.Branch("deltaNLL", deltaNLL, "deltaNLL/D")
  t.Branch("nll0", nll0, "nll0/D")
  t.Branch("nll", nll, "nll/D")
  for i, x, y in zip(range(-100, line.GetN()-100), line.GetX(), line.GetY()):
    CMS_zz4l_fai1[0] = i/100.
    deltaNLL[0] = g.Interpolate(x, y)
    t.Fill()
  f.Write()
  f.Close()
  return filename

if __name__ == "__main__":
  plotcontactterms("fL1fL1Zg_DL1_DL1Zgint_firsttest")
  plotcontactterms("fL1fL1Zg_DeR_DeLeR_firsttest")
  plotcontactterms("fL1fL1Zg_DeR_DeLint_firsttest")
  plotcontactterms("fL1fL1Zg_m1_m2_firsttest")
  plotcontactterms("fL1fL1Zg_m1_phi_firsttest")
  plotcontactterms("fL1fL1Zg_m2_phi_firsttest")
