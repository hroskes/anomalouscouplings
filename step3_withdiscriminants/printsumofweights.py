import collections
from helperstuff.enums import hypotheses, flavors
from helperstuff.samples import Sample
from math import sqrt
import os
import ROOT

samples = [Sample("qqZZ")]

sum = collections.Counter()

for sample in samples:
  f = ROOT.TFile(sample.withdiscriminantsfile())
  t = f.candTree
  for event in t:
    if t.ZZMass > 70:
      sum[sample,t.Z1Flav*t.Z2Flav] += getattr(t, sample.weightname())

for (k1, k2), v in sum.iteritems():
  print k1, k2, v
