import collections
from helperstuff.enums import hypotheses
from helperstuff.samples import Sample
from math import sqrt
import os
import ROOT

samples = [Sample("ggH", hypothesis) for hypothesis in hypotheses]

sumofweights = collections.Counter()
sumw2 = collections.Counter()

for fromsample in samples:
  t = ROOT.TChain("candTree")
  t.Add(fromsample.withdiscriminantsfile())

  length = t.GetEntries()
  for i, entry in enumerate(t, start=1):
    for tosample in samples:
      wt = getattr(t, tosample.weightname())
      sumofweights[fromsample,tosample] += wt
      sumw2[fromsample,tosample] += wt**2
    if i%10000 == 0 or i == length:
      print i, "/", length

fmt = "{:30} "*(len(samples)+2)
print fmt.format("", "rwt from", *samples)
print fmt.format("rwt to", *[""]*(len(samples)+1))
for tosample in samples:
  print fmt.format(tosample, "", *("{:.2e} +/- {:.2e}".format(sumofweights[fromsample,tosample], sqrt(sumw2[fromsample,tosample])) for fromsample in samples))
