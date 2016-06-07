import os
import ROOT
import subprocess

f = ROOT.TFile("templates_2e2mu.root")
d = f.controlPlots
if not os.path.isdir("controlplots"):
  os.mkdir("controlplots")

exts = "png", "eps", "root", "pdf"

for key in d.GetListOfKeys():
  for ext in exts:
    key.ReadObj().SaveAs("controlplots/{}.{}".format(key.ReadObj().GetName(), ext))

subprocess.call("rsync -az controlplots/ hroskes@lxplus.cern.ch:/afs/cern.ch/user/h/hroskes/www/TEST/controlplots/", shell=True)
