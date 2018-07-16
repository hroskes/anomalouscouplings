#!/usr/bin/env python

import random

def count(productionmode, rwtfrom, rwtto, category, channel, production, analysis):
  with TFile(Sample(productionmode, rwtfrom, production).withdiscriminantsfile()) if rwtfrom != "all" else DummyContextManager() as f:
    if rwtfrom == "all":
      tree = ROOT.TChain("candTree")
      effectiveentries = ROOT.TChain("effectiveentries")
      for rwtfrom in "0+", "0-", "a2", "L1", "fa3prod0.5", "fa2prod0.5", "fL1prod0.5":
        for _ in tree, effectiveentries:
          _.Add(Sample(productionmode, rwtfrom, production).withdiscriminantsfile())
    else:
      tree = f.candTree
      effectiveentries = f.effectiveentries
    t = Template(productionmode, rwtto, category, channel, production, analysis)
    hname = "".join(random.choice("abcdefghijklmnopqrstuvwxyz") for _ in xrange(25))
    tree.Draw("1>>"+hname, "({})*({})".format(t.selection, t.weightname()))
    return getattr(ROOT, hname).Integral() / sum(getattr(entry, t.weightname()) for entry in effectiveentries)

def allcats(**kwargs):
  counts = {_: count(category=_, **kwargs) for _ in categories}
  return counts

def final(**kwargs):
  print "{:12} {:7} {:12} {:12} {:12}".format("", "", *categories)
  productionmode = kwargs["productionmode"]
  rwtto = kwargs["rwtto"]
  production = kwargs["production"]
  for rwtfrom in "0+", "0-", "a2", "L1", "fa3prod0.5", "fa2prod0.5", "fL1prod0.5", "all":
      counts = allcats(rwtfrom=rwtfrom, **kwargs)
      fmt = "{:12} {:7.2g} {:>12.3g} {:>12.3g} {:>12.3g}"
      print fmt.format(
        rwtfrom,
        Sample.effectiveentries(
          Sample(productionmode, rwtfrom, production),
          ReweightingSample(rwtto, productionmode)
        ) if rwtfrom != "all" else float("nan"),
        *list(counts[category] for category in categories)
      )

import argparse
p = argparse.ArgumentParser()
p.add_argument("-a", "--analysis", choices="fa3 fa2 fL1 fL1Zg".split(), required=True)
p.add_argument("-c", "--channel", choices="2e2mu 4e 4mu".split(), required=True)
p.add_argument("-p", "--production", choices="180530 180531".split(), required=True)
p.add_argument("-r", "--rwtto", choices="0+ 0- a2 L1 L1Zg".split(), required=True)
args = p.parse_args()

from helperstuff.samples import Sample, ReweightingSample
from helperstuff.templates import Template
from helperstuff.utilities import TFile, DummyContextManager
import ROOT, rootoverloads
from helperstuff.enums import categories
c = ROOT.TCanvas()


final(productionmode="VBF", **args.__dict__)
