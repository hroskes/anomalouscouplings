#!/usr/bin/env python

from collections import Counter
from helperstuff import config
from helperstuff.enums import prodonlyhypotheses
from helperstuff.samples import Sample
import ROOT

def particletype(id):
    if abs(id) in {11, 13}: return "l"
    if abs(id) == 15: return "t"
    if id == 0: return "x"
    assert False

def particletypes(iterable):
    return Counter(particletype(id) for id in iterable)

def eventtype(leps, assocleps):
    leps = particletypes(leps)
    assocleps = particletypes(assocleps)
    if leps == Counter("llll"):
        ZZ = "4l"
    if leps == Counter("lltt"):
        ZZ = "2l2t"
    if leps == Counter("llxx"):
        ZZ = "2l2x"
    if leps == Counter("tttt"):
        ZZ = "4t"
    if leps == Counter("ttxx"):
        ZZ = "2t2x"
    if leps == Counter("xxxx"):
        ZZ = "4x"

    if assocleps == Counter("ll"):
        Z = "2l"
    if assocleps == Counter("tt"):
        Z = "2t"
    if assocleps == Counter("xx"):
        Z = "2x"

    return "ZZ{}Z{}".format(ZZ, Z)

def ZH_misreconstructed(hypothesis, production):
    s = Sample("ZH", hypothesis, production)
    t = ROOT.TChain("candTree", "candTree")
    t.Add(s.withdiscriminantsfile())
    c = Counter()
    length = t.GetEntries()
    for i, event in enumerate(t, start=1):
        if config.m4lmin <= t.ZZMass <= config.m4lmax:
            leps = (getattr(t, "GenLep{}Id".format(i)) for i in (1, 2, 3, 4))
            assocleps = (getattr(t, "GenAssocLep{}Id".format(i)) for i in (1, 2))
            c[eventtype(leps, assocleps)] += 1
        if i % 10000 == 0 or i == length: print i, "/", length
    return c

if __name__ == "__main__":
    assert len(config.productionsforcombine) == 1
    for h in prodonlyhypotheses:
        c = ZH_misreconstructed(h, config.productionsforcombine[0])
        print "{}:".format(h)
        for k, v in c.iteritems():
            print "    {} {}".format(k.replace("ZZ4", "  4").replace("ZZ2", "2").replace("Z2", " 2"), v)
