#!/usr/bin/env python
import os
import ROOT

delete = []

for filename in os.listdir("."):
  if not filename.endswith(".root"): continue
  f = ROOT.TFile(filename)
  if len(f.GetListOfKeys()) == 0 and KeepWhileOpenFile(filename+".tmp").wouldbevalid and KeepWhileOpenFile(filename.replace("_firststep", "")+".tmp").wouldbevalid:
    delete.append(filename)

if delete:
  print "rm", "\\\n   ".join(delete)
