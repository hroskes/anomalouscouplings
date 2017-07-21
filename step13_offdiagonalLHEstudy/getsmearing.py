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

hlhegenelectronpt = ROOT.TH1F("hlhegenelectronpt", "", 100, -10, 10)
hlhegenelectroneta = ROOT.TH1F("hlhegenelectroneta", "", 100, -.2, .2)
hlhegenelectronphi = ROOT.TH1F("hlhegenelectronphi", "", 100, -.2, .2)

hgenrecoelectronpt = ROOT.TH1F("hgenrecoelectronpt", "", 100, -10, 10)
hgenrecoelectroneta = ROOT.TH1F("hgenrecoelectroneta", "", 100, -.2, .2)
hgenrecoelectronphi = ROOT.TH1F("hgenrecoelectronphi", "", 100, -.2, .2)

hlherecoelectronpt = ROOT.TH1F("hlherecoelectronpt", "", 100, -10, 10)
hlherecoelectroneta = ROOT.TH1F("hlherecoelectroneta", "", 100, -.2, .2)
hlherecoelectronphi = ROOT.TH1F("hlherecoelectronphi", "", 100, -.2, .2)

hlhegenmuonpt = ROOT.TH1F("hlhegenmuonpt", "", 100, -10, 10)
hlhegenmuoneta = ROOT.TH1F("hlhegenmuoneta", "", 100, -.2, .2)
hlhegenmuonphi = ROOT.TH1F("hlhegenmuonphi", "", 100, -.2, .2)

hgenrecomuonpt = ROOT.TH1F("hgenrecomuonpt", "", 100, -10, 10)
hgenrecomuoneta = ROOT.TH1F("hgenrecomuoneta", "", 100, -.2, .2)
hgenrecomuonphi = ROOT.TH1F("hgenrecomuonphi", "", 100, -.2, .2)

hlherecomuonpt = ROOT.TH1F("hlherecomuonpt", "", 100, -10, 10)
hlherecomuoneta = ROOT.TH1F("hlherecomuoneta", "", 100, -.2, .2)
hlherecomuonphi = ROOT.TH1F("hlherecomuonphi", "", 100, -.2, .2)

hists = [
         hlherecojetpt, hlherecojeteta, hlherecojetphi,
         hlhegenelectronpt, hlhegenelectroneta, hlhegenelectronphi,
         hgenrecoelectronpt, hgenrecoelectroneta, hgenrecoelectronphi,
         hlherecoelectronpt, hlherecoelectroneta, hlherecoelectronphi,
         hlhegenmuonpt, hlhegenmuoneta, hlhegenmuonphi,
         hgenrecomuonpt, hgenrecomuoneta, hgenrecomuonphi,
         hlherecomuonpt, hlherecomuoneta, hlherecomuonphi,
        ]

length = t.GetEntries()
for i, entry in enumerate(t, start=1):
    jets = []
    LHEjets = []
    electrons = []
    genelectrons = []
    LHEelectrons = []
    muons = []
    genmuons = []
    LHEmuons = []

    for pt, eta, phi, m, id in zip(t.LHEDaughterPt, t.LHEDaughterEta, t.LHEDaughterPhi, t.LHEDaughterMass, t.LHEDaughterId):
        if abs(id) == 11:
            LHEelectrons.append(tlvfromptetaphim(pt, eta, phi, m))
        elif abs(id) == 13:
            LHEmuons.append(tlvfromptetaphim(pt, eta, phi, m))
    for pt, eta, phi, m, id in zip(t.LHEAssociatedParticlePt, t.LHEAssociatedParticleEta, t.LHEAssociatedParticlePhi, t.LHEAssociatedParticleMass, t.LHEAssociatedParticleId):
        if 1 <= abs(id) <= 6 or id == 21:
            LHEjets.append(tlvfromptetaphim(pt, eta, phi, m))
    for pt, eta, phi, id in zip(*[[getattr(t, "GenLep{}{}".format(j, var)) for j in range(1, 5)] for var in ("Pt", "Eta", "Phi", "Id")]):
        m = 0
        if abs(id) == 11:
            genelectrons.append(tlvfromptetaphim(pt, eta, phi, m))
        elif abs(id) == 13:
            genmuons.append(tlvfromptetaphim(pt, eta, phi, m))
    for pt, eta, phi, id in zip(t.LepPt, t.LepEta, t.LepPhi, t.LepLepId):
        if abs(id) == 11:
            electrons.append(tlvfromptetaphim(pt, eta, phi, m))
        elif abs(id) == 13:
            muons.append(tlvfromptetaphim(pt, eta, phi, m))
    for pt, eta, phi, mass in zip(t.JetPt, t.JetEta, t.JetPhi, t.JetMass):
        jets.append(tlvfromptetaphim(pt, eta, phi, 0))

    for lhejet in LHEjets:
        if not jets: continue
        recojet = min(jets, key=lambda jet: jet.DeltaR(lhejet))
        if lhejet != min(LHEjets, key=lambda jet: jet.DeltaR(recojet)): continue
        hlherecojetpt.Fill(recojet.Pt() - lhejet.Pt())
        hlherecojeteta.Fill(recojet.Eta() - lhejet.Eta())
        hlherecojetphi.Fill(recojet.Phi() - lhejet.Phi())

    for lheelectron in LHEelectrons:
        recoelectron = min(electrons, key=lambda electron: electron.DeltaR(lheelectron))
        if lheelectron != min(LHEelectrons, key=lambda electron: electron.DeltaR(recoelectron)): continue
        hlherecoelectronpt.Fill(recoelectron.Pt() - lheelectron.Pt())
        hlherecoelectroneta.Fill(recoelectron.Eta() - lheelectron.Eta())
        hlherecoelectronphi.Fill(recoelectron.Phi() - lheelectron.Phi())

    for genelectron in genelectrons:
        recoelectron = min(electrons, key=lambda electron: electron.DeltaR(genelectron))
        if genelectron != min(genelectrons, key=lambda electron: electron.DeltaR(recoelectron)): continue
        hgenrecoelectronpt.Fill(recoelectron.Pt() - genelectron.Pt())
        hgenrecoelectroneta.Fill(recoelectron.Eta() - genelectron.Eta())
        hgenrecoelectronphi.Fill(recoelectron.Phi() - genelectron.Phi())

    for lheelectron in LHEelectrons:
        genelectron = min(genelectrons, key=lambda electron: electron.DeltaR(lheelectron))
        if lheelectron != min(LHEelectrons, key=lambda electron: electron.DeltaR(genelectron)): continue
        hlhegenelectronpt.Fill(genelectron.Pt() - lheelectron.Pt())
        hlhegenelectroneta.Fill(genelectron.Eta() - lheelectron.Eta())
        hlhegenelectronphi.Fill(genelectron.Phi() - lheelectron.Phi())

    for lhemuon in LHEmuons:
        recomuon = min(muons, key=lambda muon: muon.DeltaR(lhemuon))
        if lhemuon != min(LHEmuons, key=lambda muon: muon.DeltaR(recomuon)): continue
        hlherecomuonpt.Fill(recomuon.Pt() - lhemuon.Pt())
        hlherecomuoneta.Fill(recomuon.Eta() - lhemuon.Eta())
        hlherecomuonphi.Fill(recomuon.Phi() - lhemuon.Phi())

    for genmuon in genmuons:
        recomuon = min(muons, key=lambda muon: muon.DeltaR(genmuon))
        if genmuon != min(genmuons, key=lambda muon: muon.DeltaR(recomuon)): continue
        hgenrecomuonpt.Fill(recomuon.Pt() - genmuon.Pt())
        hgenrecomuoneta.Fill(recomuon.Eta() - genmuon.Eta())
        hgenrecomuonphi.Fill(recomuon.Phi() - genmuon.Phi())

    for lhemuon in LHEmuons:
        genmuon = min(genmuons, key=lambda muon: muon.DeltaR(lhemuon))
        if lhemuon != min(LHEmuons, key=lambda muon: muon.DeltaR(genmuon)): continue
        hlhegenmuonpt.Fill(genmuon.Pt() - lhemuon.Pt())
        hlhegenmuoneta.Fill(genmuon.Eta() - lhemuon.Eta())
        hlhegenmuonphi.Fill(genmuon.Phi() - lhemuon.Phi())

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
