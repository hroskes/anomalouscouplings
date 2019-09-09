#!/usr/bin/env python
import argparse

p = argparse.ArgumentParser()
p.add_argument("directory", nargs="*", type=lambda x: x.rstrip("/"))
args = p.parse_args()

import os

from helperstuff.utilities import cd, KeepWhileOpenFile, TFile

delete = []

with cd(os.path.dirname(__file__)):
  for directory in os.listdir("."):
    if not os.path.isdir(directory): continue
    if args.directory and directory not in args.directory: continue
    for filename in os.listdir(directory):
      filename = os.path.join(directory, filename)
      if not filename.endswith(".root"): continue
      if not KeepWhileOpenFile(filename+".tmp").wouldbevalid: continue
      try:
        with TFile(filename) as f:
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
