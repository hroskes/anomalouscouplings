#!/usr/bin/env python
from collections import Counter
from helperstuff import config
from helperstuff.combinehelpers import getrate
from helperstuff.enums import categories, Channel, channels, ProductionMode
from helperstuff.templates import DataTree
from helperstuff.utilities import tfiles
import os
from projections import Projections
import ROOT
import sys

if __name__ == "__main__":
  totalexp = totalobs = 0
  exp = Counter()
  for production in config.productionsforcombine:
    print production
    productionmodes = [ProductionMode(p) for p in ("ggH", "VBF", "ZH", "WH", "ttH", "qqZZ", "ggZZ", "VBF bkg", "ZX")]
    for category in categories:
      for flavor in channels:
        print flavor
        print
        for p in productionmodes:
          exp[p,flavor,category] += getrate(p, flavor, "forexpectedscan", production, category, *sys.argv[1:])
          totalexp += exp[p,flavor,category]
        print " ".join("{:.2f}" for p in productionmodes).format(*(exp[p,flavor,category] for p in productionmodes))
        print "Total signal:   {:.2f}".format(sum(exp[p,flavor,category] for p in productionmodes if not p.isbkg))
        print "Total bkg:      {:.2f}".format(sum(exp[p,flavor,category] for p in productionmodes if p.isbkg))
        print "Total expected: {:.2f}".format(sum(exp[p,flavor,category] for p in productionmodes))
        print "Observed:", tfiles[DataTree(production, flavor, category).treefile].candTree.GetEntries()
        totalobs += tfiles[DataTree(production, flavor, category).treefile].candTree.GetEntries()
        print
    print

  print
  print
  print "Table:"
  print
  flavors = [Channel(f) for f in "4e", "4mu", "2e2mu"]
  fmt = "{:<10}"*5
  print fmt.format(*([""]+flavors+["total"]))
  for c in categories:
    totalsig = Counter()
    totalbkg = Counter()
    print c
    for p in productionmodes:
      print fmt.format(*([p] + ["{:.2f}".format(exp[p,f,c]) for f in flavors] + ["{:.2f}".format(sum(exp[p,f,c] for f in flavors))]))
      for f in flavors:
        if p.isbkg:
          totalbkg[f] += exp[p,f,c]
        else:
          totalsig[f] += exp[p,f,c]
    print
    print fmt.format(*(["bkg"] + ["{:.2f}".format(totalbkg[f]) for f in flavors] + ["{:.2f}".format(sum(totalbkg[f] for f in flavors))]))
    print fmt.format(*(["signal"] + ["{:.2f}".format(totalsig[f]) for f in flavors] + ["{:.2f}".format(sum(totalsig[f] for f in flavors))]))
    print fmt.format(*(["observed"] + [sum(tfiles[DataTree(production, flavor, category).treefile].candTree.GetEntries() for production in config.productionsforcombine) for flavor in flavors] + [""]))
    print
    print
  print "Total:"
  print "Expected:", totalexp
  print "Observed:", totalobs
