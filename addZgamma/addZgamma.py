#!/usr/bin/env python

import array
import os
import sys

import ROOT

from ZZMatrixElement.PythonWrapper.mela import Mela, SimpleParticle_t, SimpleParticleCollection_t, TVar

mela = Mela(125, 13)

def leptons(t):
    leptons = SimpleParticleCollection_t()
    for id, pt, eta, phi in zip(t.LepLepId, t.LepPt, t.LepEta, t.LepPhi):
        m = 0
        tlv = ROOT.TLorentzVector()
        tlv.SetPtEtaPhiM(pt, eta, phi, m)
        leptons.push_back(SimpleParticle_t(id, tlv))
    return leptons

def addL1Zg(t, doreweighting=False):
    newtree = t.CloneTree(0)
    SM = array.array('f', [0])
    L1Zg = array.array('f', [0])
    L1Zgint = array.array('f', [0])
    a3test = array.array('f', [0])
    newtree.Branch("p0plus_forL1Zg_VAJHU", a3test, "p0plus_forL1Zg_VAJHU/F")
    newtree.Branch("p0minus_forL1Zg_VAJHU", SM, "p0minus_forL1Zg_VAJHU/F")
    newtree.Branch("p0_ghzgs1prime2_VAJHU", L1Zg, "p0_ghzgs1prime2_VAJHU/F")
    newtree.Branch("pg1ghzgs1prime2_VAJHU", L1Zgint, "pg1ghzgs1prime2_VAJHU/F")

    for entry in t:
        mela.setInputEvent(leptons(t), 0, 0, False)

        mela.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZINDEPENDENT)
        mela.ghz1 = 1
        SM[0] = mela.computeP(True)

        mela.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZINDEPENDENT)
        mela.ghz4 = 1
        a3test[0] = mela.computeP(True)

        mela.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZINDEPENDENT)
        mela.ghzgs1_prime2 = 1
        L1Zg[0] = mela.computeP(True)

        mela.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.ZZINDEPENDENT)
        mela.ghz1 = mela.ghzgs1_prime2 = 1
        L1Zgint[0] = mela.computeP(True)

        mela.resetInputEvent()
        newtree.Fill()

    return newtree

if __name__ == "__main__":
    try:
        infile = sys.argv[1]
        outfile = sys.argv[2]
        if not (infile.endswith(".root") and os.path.exists(infile)):
            raise ValueError("infile should end with .root and exist!")
        if not (outfile.endswith(".root") and not os.path.exists(outfile)):
            raise ValueError("outfile should end with .root and not exist!")
    except:
        print "python", sys.argv[0], "infile.root outfile.root"
        raise

    f = ROOT.TFile(infile)
    t = f.Get("ZZTree/candTree")
    try:
        newf = ROOT.TFile(outfile, "RECREATE")
        newdir = newf.mkdir("ZZTree")
        newdir.cd()
        newt = addL1Zg(t)
        newt.Write()
    except:
        try:
            raise
        finally:
            try:
                os.remove(outfile)
            except:
                pass
