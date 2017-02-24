#!/usr/bin/env python
import os
import ROOT

delete = []

for filename in os.listdir("."):
  if not filename.endswith(".root"): continue
  if os.path.exists(filename+".tmp"): continue
  f = ROOT.TFile(filename)
  try:
    if f.candTree.GetEntries() == 0:
      delete.append(filename)
  except AttributeError:
      delete.append(filename)

  if filename in delete and os.path.exists(filename+".tmp"): delete.remove(filename)

if delete:
    print "rm \t\t\t\\"
    print "\t\t\t\\\n".join(delete)
    print
