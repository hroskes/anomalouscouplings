#!/usr/bin/env python
import os

import ROOT

from helperstuff import config

from helperstuff.samples import Sample
from helperstuff.utilities import cache, mkdir_p, tfiles, tlvfromptetaphim

TF1 = cache(ROOT.TF1)

s = Sample("VBF", "0+", config.productionforcombine)

f = tfiles[s.withdiscriminantsfile()]
t = f.candTree

hlherecojetpt = ROOT.TH1F("hlherecojetpt", "", 100, -100, 100)
hlherecojeteta = ROOT.TH1F("hlherecojeteta", "", 100, -1, 1)
hlherecojetphi = ROOT.TH1F("hlherecojetphi", "", 100, -1, 1)

hlhegenleptonpt = ROOT.TH1F("hlhegenleptonpt", "", 100, -10, 10)
hlhegenleptoneta = ROOT.TH1F("hlhegenleptoneta", "", 100, -.2, .2)
hlhegenleptonphi = ROOT.TH1F("hlhegenleptonphi", "", 100, -.2, .2)

hgenrecoleptonpt = ROOT.TH1F("hgenrecoleptonpt", "", 100, -10, 10)
hgenrecoleptoneta = ROOT.TH1F("hgenrecoleptoneta", "", 100, -.2, .2)
hgenrecoleptonphi = ROOT.TH1F("hgenrecoleptonphi", "", 100, -.2, .2)

hlherecoleptonpt = ROOT.TH1F("hlherecoleptonpt", "", 100, -10, 10)
hlherecoleptoneta = ROOT.TH1F("hlherecoleptoneta", "", 100, -.2, .2)
hlherecoleptonphi = ROOT.TH1F("hlherecoleptonphi", "", 100, -.2, .2)

hists = [hlherecojetpt, hlherecojeteta, hlherecojetphi, hlhegenleptonpt, hlhegenleptoneta, hlhegenleptonphi, hgenrecoleptonpt, hgenrecoleptoneta, hgenrecoleptonphi, hlherecoleptonpt, hlherecoleptoneta, hlherecoleptonphi]

length = t.GetEntries()
for i, entry in enumerate(t, start=1):
    jets = []
    LHEjets = []
    leptons = []
    genleptons = []
    LHEleptons = []

    for pt, eta, phi, m, id in zip(t.LHEDaughterPt, t.LHEDaughterEta, t.LHEDaughterPhi, t.LHEDaughterMass, t.LHEDaughterId):
        LHEleptons.append(tlvfromptetaphim(pt, eta, phi, m))
    for pt, eta, phi, m, id in zip(t.LHEAssociatedParticlePt, t.LHEAssociatedParticleEta, t.LHEAssociatedParticlePhi, t.LHEAssociatedParticleMass, t.LHEAssociatedParticleId):
        if 1 <= abs(id) <= 6 or id == 21:
            LHEjets.append(tlvfromptetaphim(pt, eta, phi, m))
    for pt, eta, phi, id in zip(*[[getattr(t, "GenLep{}{}".format(j, var)) for j in range(1, 5)] for var in ("Pt", "Eta", "Phi", "Id")]):
        m = 0
        genleptons.append(tlvfromptetaphim(pt, eta, phi, m))
    for pt, eta, phi, id in zip(t.LepPt, t.LepEta, t.LepPhi, t.LepLepId):
        leptons.append(tlvfromptetaphim(pt, eta, phi, 0))
    for pt, eta, phi, mass in zip(t.JetPt, t.JetEta, t.JetPhi, t.JetMass):
        jets.append(tlvfromptetaphim(pt, eta, phi, 0))

    for lhejet in LHEjets:
        if not jets: continue
        recojet = min(jets, key=lambda jet: jet.DeltaR(lhejet))
        if lhejet != min(LHEjets, key=lambda jet: jet.DeltaR(recojet)): continue
        hlherecojetpt.Fill(recojet.Pt() - lhejet.Pt())
        hlherecojeteta.Fill(recojet.Eta() - lhejet.Eta())
        hlherecojetphi.Fill(recojet.Phi() - lhejet.Phi())

    for lhelepton in LHEleptons:
        recolepton = min(leptons, key=lambda lepton: lepton.DeltaR(lhelepton))
        if lhelepton != min(LHEleptons, key=lambda lepton: lepton.DeltaR(recolepton)): continue
        hlherecoleptonpt.Fill(recolepton.Pt() - lhelepton.Pt())
        hlherecoleptoneta.Fill(recolepton.Eta() - lhelepton.Eta())
        hlherecoleptonphi.Fill(recolepton.Phi() - lhelepton.Phi())

    for genlepton in genleptons:
        recolepton = min(leptons, key=lambda lepton: lepton.DeltaR(genlepton))
        if genlepton != min(genleptons, key=lambda lepton: lepton.DeltaR(recolepton)): continue
        hgenrecoleptonpt.Fill(recolepton.Pt() - genlepton.Pt())
        hgenrecoleptoneta.Fill(recolepton.Eta() - genlepton.Eta())
        hgenrecoleptonphi.Fill(recolepton.Phi() - genlepton.Phi())

    for lhelepton in LHEleptons:
        genlepton = min(genleptons, key=lambda lepton: lepton.DeltaR(lhelepton))
        if lhelepton != min(LHEleptons, key=lambda lepton: lepton.DeltaR(genlepton)): continue
        hlhegenleptonpt.Fill(genlepton.Pt() - lhelepton.Pt())
        hlhegenleptoneta.Fill(genlepton.Eta() - lhelepton.Eta())
        hlhegenleptonphi.Fill(genlepton.Phi() - lhelepton.Phi())

    if i % 1000 == 0 or i == length:
        print i, "/", length

c = ROOT.TCanvas()

saveasdir = os.path.join(config.plotsbasedir, "offdiagonalLHEstudy", "resolution")
mkdir_p(saveasdir)
for h in hists:
    for ext in "png eps root pdf".split():
        h.Draw()
        f = TF1("f"+h.GetName(), "gaus(0)", h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
        f.SetParameters(h.GetEntries(), h.GetMean(), h.GetRMS())
        h.Fit(f)
        c.SaveAs(os.path.join(saveasdir, h.GetName()+"."+ext))
