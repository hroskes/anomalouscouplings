import ROOT

from ZZMatrixElement.PythonWrapper.mela import Mela, SimpleParticle_t, SimpleParticleCollection_t, TVar

mela = Mela()

def leptons(t):
    leptons = SimpleParticleCollection_t()
    for i in range(1, 5):
        id, pt, eta, phi = (getattr(t, "GenLep{}{}".format(i, thing)) for thing in ("Id", "Pt", "Eta", "Phi"))
        m = 0
        tlv = ROOT.TLorentzVector()
        tlv.SetPtEtaPhiM(pt, eta, phi, m)
        leptons.push_back(SimpleParticle_t(id, tlv)
    return leptons

def addL1Zg(t, doreweighting=False):
    newtree = t.Clone(0)
    L1Zg = array.array('f', [0])
    L1Zgint = array.array('f', [0])
    newtree.Branch("p0_ghzgs1prime2_VAJHU", L1Zg, "p0_ghzgs1prime2_VAJHU/F")
    newtree.Branch("pg1ghzgs1prime2_VAJHU", L1Zgint, "pg1ghzgs1prime2_VAJHU/F")

    for entry in t:
        mela.setInputEvent(leptons(t), 0, 0, False)
        mela.ghzgs1prime2 = mela.ghg2 = 1
        L1Zg[0] = mela.computeP(True)
        mela.ghz1 = mela.ghg2 = mela.ghzgs1prime2 = 1
        L1Zgint[0] = mela.computeP(True)
        mela.resetInputEvent()
        newtree.Fill()

    return newtree
