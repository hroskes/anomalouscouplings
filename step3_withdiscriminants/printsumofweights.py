import collections
from helperstuff import config
from helperstuff.enums import prodonlyhypotheses, flavors, Channel, Production
from helperstuff.samples import ReweightingSample, Sample
from math import sqrt
import os
import ROOT

samples = [ReweightingSample("VBF", h) for h in prodonlyhypotheses]
production = Production("160729")

flavordict = {13**4: Channel("4mu"), 11**4: Channel("4e"), 11**2*13**2: Channel("2e2mu")}

#https://docs.python.org/2/library/collections.html
class OrderedCounter(collections.Counter, collections.OrderedDict):
     'Counter that remembers the order elements are first encountered'

     def __repr__(self):
         return '%s(%r)' % (self.__class__.__name__, OrderedDict(self))

     def __reduce__(self):
         return self.__class__, (OrderedDict(self),)

sum = OrderedCounter()
for reweightsample in samples:
  for sample in samples:
    sum[sample,reweightsample] = 0  #to lock in the order

for sample in samples:
  f = ROOT.TFile(Sample(sample, production).withdiscriminantsfile())
  t = f.candTree
  length = t.GetEntries()
  for i, event in enumerate(t, start=1):
    if config.m4lmin < t.ZZMass < config.m4lmax:
      for reweightsample in samples:#sample.reweightingsamples():
        #sum[sample,reweightsample,flavordict[t.Z1Flav*t.Z2Flav]] += getattr(t, reweightsample.weightname()) * production.dataluminosity
        sum[sample,reweightsample] += getattr(t, reweightsample.weightname()) * production.dataluminosity
    if i % 10000 == 0:
      print i, "/", length

for (k1, k2), v in sum.iteritems():
  print "{:<20} --> {:<20} {}".format(k1, k2, v)
