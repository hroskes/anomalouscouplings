#!/usr/bin/env python
import os
import ROOT

delete = []

for filename in os.listdir("."):
  if not filename.endswith(".root"): continue
  f = ROOT.TFile(filename)
  if len(f.GetListOfKeys()) == 0 and not os.path.exists(filename+".tmp"):
    delete.append(filename)

print "rm", "\\\n   ".join(delete)
