#!/usr/bin/env python

assert __name__ == "__main__"

import argparse
import os

p = argparse.ArgumentParser()
p.add_argument("folder", nargs="*", default=os.listdir(os.path.dirname(__file__)))
args = p.parse_args()

import ROOT
from helperstuff.utilities import cd, KeepWhileOpenFile

delete = []

with cd(os.path.dirname(__file__)):
  for directory in args.folder:
    if not os.path.isdir(directory): continue
    for filename in os.listdir(directory):
      filename = os.path.join(directory, filename)
      if not filename.endswith(".root"): continue
      if not KeepWhileOpenFile(filename+".tmp").wouldbevalid: continue
      if not KeepWhileOpenFile(filename.replace("_firststep", "")+".tmp").wouldbevalid: continue
      f = ROOT.TFile(filename)
      try:
        if len(f.GetListOfKeys()) == 0:
          delete.append(filename)
      except AttributeError:
        delete.append(filename)

for filename in delete[:]:
  if not KeepWhileOpenFile(filename+".tmp").wouldbevalid: delete.remove(filename)

if delete:
  print "rm", "\\\n   ".join(delete)
