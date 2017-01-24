#!/usr/bin/env python

import ROOT

from helperstuff.CJLSTscripts import *
from helperstuff.utilities import *

f = ROOT.TFile("VBF0+170119.root")
t = f.candTree

h = ROOT.TH1F("h", "h", 100, 40, 220)
hdelta = ROOT.TH1F("h", "h", 100, -200, 200)


length = t.GetEntries()
for i, entry in enumerate(t, start=1):
#    if t.category_0P != VHHadrTaggedIchep16: continue
    jets = [tlvfromptetaphim(pt, eta, phi, m) for pt, eta, phi, m in zip(t.LHEAssociatedParticlePt, t.LHEAssociatedParticleEta, t.LHEAssociatedParticlePhi, t.LHEAssociatedParticleMass)]
    if len(jets) < 2: continue
    mjj = (jets[0]+jets[1]).M()
    h.Fill(mjj)
    hdelta.Fill(t.DiJetMass-mjj)
    if i % 1000 == 0 or i == length:
       print i, "/", length

c = ROOT.TCanvas()
h.Draw()
c.SaveAs("~/www/TEST/mjj_lhe_all.png")
hdelta.Draw()
c.SaveAs("~/www/TEST/deltamjj_all.png")

