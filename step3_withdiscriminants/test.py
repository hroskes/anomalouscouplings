#!/usr/bin/env python

import glob, ROOT, rootoverloads

for fi in glob.iglob("190821_201*/*.root"):
  f = ROOT.TFile(fi)
  t = f.candTree
  t.SetBranchStatus("*", 0)
  t.SetBranchStatus("D_2jet_0plus", 1)
  t.GetEntry(t.GetEntries() - 1)
  t.D_2jet_0plus
