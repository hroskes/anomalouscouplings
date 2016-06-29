import collections
from helperstuff.enums import hypotheses, flavors, Channel
from helperstuff.samples import Sample
from math import sqrt
import os
import ROOT

samples = [Sample("ggH", "fa30.5")]

flavordict = {13**4: Channel("4mu"), 11**4: Channel("4e"), 11**2*13**2: Channel("2e2mu")}

sum = collections.Counter()

for sample in samples:
  f = ROOT.TFile(sample.withdiscriminantsfile())
  t = f.candTree
  length = t.GetEntries()
  for i, event in enumerate(t, start=1):
    if 105 < t.ZZMass < 140:
      for reweightsample in sample.reweightingsamples():
        sum[sample,reweightsample,flavordict[t.Z1Flav*t.Z2Flav]] += getattr(t, reweightsample.weightname())
    if i % 10000 == 0:
      print i, "/", length

for (k1, k2, k3), v in sum.iteritems():
  print k1, "-->", k2, k3, v
