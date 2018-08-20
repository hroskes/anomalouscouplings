#!/usr/bin/env python
import os

import ROOT

from helperstuff.utilities import KeepWhileOpenFile

delete = []

for filename in os.listdir("."):
  if not filename.endswith(".root"): continue
  if not KeepWhileOpenFile(filename+".tmp").wouldbevalid: continue
  f = ROOT.TFile(filename)
  try:
    if f.candTree.GetEntries() == 0:
      delete.append(filename)
  except AttributeError:
      delete.append(filename)

  if filename in delete and KeepWhileOpenFile(filename+".tmp").wouldbevalid: delete.remove(filename)

if delete:
    print "rm \t\t\t\\"
    print "\t\t\t\\\n".join(delete)
    print
