import os
import ROOT

delete = []

for filename in os.listdir("."):
  if not filename.endswith(".root"): continue
  f = ROOT.TFile(filename)
  if len(f.GetListOfKeys()) == 0:
    delete.append(filename)

print "rm", " ".join(delete)
