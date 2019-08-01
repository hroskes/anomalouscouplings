#!/usr/bin/env python

import ROOT

from helperstuff import config
from helperstuff.CJLSTscripts import VHHadrTaggedMor18, VBF2jTaggedMor18, getDVBF2jetsConstant, getDVBF2jetsWP
from helperstuff.samples import Sample
from helperstuff.utilities import tfiles

assert len(config.productionsforcombine) == 1
fromsample = Sample("VBF", "0+", config.productionsforcombine[0])
tosample = fromsample
categorization = "category_0P"

t = tfiles[fromsample.withdiscriminantsfile()].candTree

weightname = tosample.MC_weight
n = n_total = 0.
length = t.GetEntries()
for i, entry in enumerate(t):
    ZZMass = t.ZZMass
    if not config.m4lmin < ZZMass < config.m4lmax:
        continue
    weight = 1
    n_total += weight
    if getattr(t, categorization) != VBF2jTaggedMor18 and t.p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal > 0 and 1/(1+getDVBF2jetsConstant(ZZMass)*t.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal/t.p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal) > getDVBF2jetsWP(ZZMass, 0) and t.nExtraLep == 0:
        n += weight
    if i % 10000 == 0 or i == length:
        print i, "/", length

print n, n_total, n/n_total
