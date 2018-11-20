#!/usr/bin/env python
import os

import ROOT

from helperstuff.utilities import cd, KeepWhileOpenFile

delete = []

with cd(os.path.dirname(__file__)):
  for directory in os.listdir("."):
    if not os.path.isdir(directory): continue
    for filename in os.listdir(directory):
      filename = os.path.join(directory, filename)
      if not filename.endswith(".root"): continue
      if not KeepWhileOpenFile(filename+".tmp").wouldbevalid: continue
      f = ROOT.TFile(filename)
      try:
        if f.candTree.GetEntries() == 0:
          delete.append(filename)
      except AttributeError:
        delete.append(filename)

for filename in delete[:]:
  if not KeepWhileOpenFile(filename+".tmp").wouldbevalid: delete.remove(filename)

if delete:
  print "rm \t\t\t\\"
  print "\t\t\t\\\n".join(delete)
  print
