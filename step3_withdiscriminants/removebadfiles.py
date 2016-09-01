import os
import ROOT

delete = []

for filename in os.listdir("."):
  if not filename.endswith(".root"): continue
  f = ROOT.TFile(filename)
  try:
    if f.candTree.GetEntries() == 0:
      delete.append(filename)
  except AttributeError:
      delete.append(filename)

print "rm \t\t\t\\"
for a in delete:
    print a, "\t\t\t\\"
print
