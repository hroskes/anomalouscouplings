#!/usr/bin/env python

import sys

import ROOT

import rootoverloads

from helperstuff import config, style
from helperstuff.samples import Sample
from helperstuff.utilities import tfiles

s = Sample(config.productionforcombine, *sys.argv[1:])
t = tfiles[s.withdiscriminantsfile()].candTree

t.SetBranchStatus("*", 0)
t.SetBranchStatus("LHE*Id", 1)
t.SetBranchStatus("ZZMass", 1)

h1 = ROOT.TH1F("h1", "", 35, 105, 140)
h2 = ROOT.TH1F("h2", "", 35, 105, 140)

nentries = 0
nnotaus = 0
nmass = 0
nmassnotaus = 0
for entry in t:
  if 105 < t.ZZMass < 140: nmass += 1
  if sum(_ in (-11, 11, -13, 13) for _ in list(t.LHEDaughterId) + list(t.LHEAssociatedParticleId)) >= 4:
    nnotaus += 1
    if 105 < t.ZZMass < 140:
      nmassnotaus += 1
    h1.Fill(t.ZZMass)
  else:
    h2.Fill(t.ZZMass)
  nentries += 1

print "Total:", nentries
print "  With taus:", nentries - nnotaus
print "  No taus:", nnotaus
print "105 - 140:", nmass
print "  With taus:", nmass - nmassnotaus
print "  No taus:", nmassnotaus

print float(nmass - nmassnotaus) / nmass

hs = ROOT.THStack("hs", "hs")
l = ROOT.TLegend(.6, .6, .9, .9)
hs.Add(h1)
hs.Add(h2)
h1.SetLineColor(2)
h2.SetLineColor(4)
l.AddEntry(h1, "4 e/#mu")
l.AddEntry(h2, "no 4 e/#mu")
l.SetBorderSize(0)
l.SetFillStyle(0)

print h1.GetEntries(), h2.GetEntries()

c = ROOT.TCanvas()
hs.Draw()
l.Draw()
c.SaveAs("~/www/TEST/test.png")
