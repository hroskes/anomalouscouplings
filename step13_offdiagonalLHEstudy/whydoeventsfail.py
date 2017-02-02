#!/usr/bin/env python

import ROOT

def GenLepPt(self, i=None):
    if i is None:
        return [self.GenLepPt(i) for i in range(1, 5)]
    return getattr(self, "GenLep{}Pt".format(i))
def GenLepEta(self, i=None):
    if i is None:
        return [self.GenLepEta(i) for i in range(1, 5)]
    return getattr(self, "GenLep{}Eta".format(i))
def GenLepId(self, i=None):
    if i is None:
        return [self.GenLepId(i) for i in range(1, 5)]
    return getattr(self, "GenLep{}Id".format(i))

ROOT.TTree.GenLepPt = GenLepPt
ROOT.TTree.GenLepId = GenLepId
ROOT.TTree.GenLepEta = GenLepEta

def whydoeseventfail(t):
    if t.isSelected: return None
    for i in range(1, 5):
        if t.GenLepPt(i)<5 and abs(t.GenLepId(i))==13: return "muon pt too low {}".format(t.GenLepPt(i))
        if t.GenLepPt(i)<7 and abs(t.GenLepId(i))==11: return "electron pt too low {}".format(t.GenLepPt(i))
        if abs(t.GenLepEta(i))>2.4 and abs(t.GenLepId(i))==13: return "muon eta too high {}".format(t.GenLepEta(i))
        if abs(t.GenLepEta(i))>2.5 and abs(t.GenLepId(i))==11: return "electron eta too high {}".format(t.GenLepEta(i))
    t.Show()
    print "id: ", t.GenLepId()
    print "pt: ", t.GenLepPt()
    print "eta:", t.GenLepEta()
    assert False

if __name__ == "__main__":
    f = ROOT.TFile("/afs/cern.ch/work/w/wqin/public/forHeshy/orig_pythia.root")
    t = f.SelectedTree
    for entry in t:
        result = whydoeseventfail(t)
        if result:
            print result
