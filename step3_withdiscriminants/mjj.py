#!/usr/bin/env python

import ROOT

from helperstuff.CJLSTscripts import *
from helperstuff.utilities import *

f = ROOT.TFile("VBF0+170119.root")
t = f.candTree

category = False

h = ROOT.TH1F("h", "h", 100, 0, 7000)
if category:
    hdelta = ROOT.TH1F("hdelta", "h", 100, -2300, 300)
else:
    hdelta = ROOT.TH1F("hdelta", "h", 100, -5000, 7000)


length = t.GetEntries()
for i, entry in enumerate(t, start=1):
    if t.category_0P != VHHadrTaggedIchep16: continue
    jets = [tlvfromptetaphim(pt, eta, phi, m) for pt, eta, phi, m in zip(t.LHEAssociatedParticlePt, t.LHEAssociatedParticleEta, t.LHEAssociatedParticlePhi, t.LHEAssociatedParticleMass)]
    if len(jets) < 2: continue
    mjj = (jets[0]+jets[1]).M()
    h.Fill(mjj)
    hdelta.Fill(t.DiJetMass-mjj)
    if i % 1000 == 0 or i == length:
       print i, "/", length

append = ""
if not category:
    append = "_all"

c = ROOT.TCanvas()
h.Draw()
c.SaveAs("~/www/TEST/mjj_lhe{}.png".format(append))
hdelta.Draw()
c.SaveAs("~/www/TEST/deltamjj{}.png".format(append))

