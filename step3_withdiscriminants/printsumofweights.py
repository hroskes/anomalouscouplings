import collections
from helperstuff.enums import hypotheses, flavors, Channel
from helperstuff.samples import Sample
from math import sqrt
import os
import ROOT

samples = [Sample("qqZZ"), Sample("ggH", "0+")]

flavordict = {13**4: Channel("4mu"), 11**4: Channel("4e"), 11**2*13**2: Channel("2e2mu")}

sum = collections.Counter()

for sample in samples:
  f = ROOT.TFile(sample.withdiscriminantsfile())
  t = f.candTree
  for event in t:
    if 105 < t.ZZMass < 140:
      sum[sample,flavordict[t.Z1Flav*t.Z2Flav]] += getattr(t, sample.weightname())

for (k1, k2), v in sum.iteritems():
  print k1, k2, v
