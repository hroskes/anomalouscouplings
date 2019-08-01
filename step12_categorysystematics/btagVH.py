#!/usr/bin/env python

import ROOT

from helperstuff import config
from helperstuff.CJLSTscripts import VHHadrTaggedAC19, VBF2jTaggedAC19
from helperstuff.samples import Sample
from helperstuff.utilities import tfiles

assert len(config.productionsforcombine) == 0
fromsample = Sample("ZH", "0+", config.productionsforcombine[0], "POWHEG")
tosample = fromsample
categorization = "category_0P_or_0M"

t = tfiles[fromsample.withdiscriminantsfile()].candTree

weightparts_total = [
  tosample.MC_weight,
  "ZZMass > {}".format(config.m4lmin),
  "ZZMass < {}".format(config.m4lmax),
]
weight_total = " && ".join("("+_+")" for _ in weightparts_total)
weightparts = weightparts_total + [
  "{} != {}".format(categorization, VHHadrTaggedAC19),
  "{} != {}".format(categorization, VBF2jTaggedAC19),
  "nExtraLep==0 && (nCleanedJetsPt30==2||nCleanedJetsPt30==3) && nCleanedJetsPt30BTagged>=2",
]
weight = " && ".join("("+_+")" for _ in weightparts)

t.Draw("1>>h", weight)
n = ROOT.h.Integral()
t.Draw("1>>h_total", weight_total)
n_total = ROOT.h_total.Integral()
print n, n_total, n/n_total
