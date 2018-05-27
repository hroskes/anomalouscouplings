#!/usr/bin/env python

assert __name__ == "__main__"

import argparse

p = argparse.ArgumentParser()
p.add_argument("filename")
args = p.parse_args()

import glob
import itertools
import os

import ROOT

f = ROOT.TFile(args.filename)

t = f.candTree

print sum(1 for a in t if t.category_0P_or_0M == 2)
print sum(1 for a in t if t.category_0P_or_0M != 2 and max(t.D_2jet_0plus, t.D_2jet_0minus) > 0.5)
print sum(1 for a in t if t.category_0P_or_0M == 4)
print sum(1 for a in t if t.category_0P_or_0M not in (2, 4) and max(t.D_HadWH_0plus, t.D_HadWH_0minus, t.D_HadZH_0plus, t.D_HadZH_0minus) > 0.5)
