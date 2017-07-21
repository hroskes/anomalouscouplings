import getpass
import os
import re
import socket

def getconfiguration(hostname, username):
  if (".cern.ch" in hostname or "lxplus" in hostname) and username == "hroskes":
    return dict(
      host = "lxplus",
      repositorydir = "/afs/cern.ch/work/h/hroskes/LHEanomalouscouplings/",
      plotsbasedir = "/afs/cern.ch/user/h/hroskes/www/LHEanomalouscouplings/",
      svndir = "/afs/cern.ch/work/h/hroskes/AN/notes",
    )

  elif ("login-node" in hostname or "compute" in hostname or "bigmem" in hostname) and username == "jroskes1@jhu.edu":
    return dict(
      host = "MARCC",
      repositorydir = "/work-zfs/lhc/heshy/LHEanomalouscouplings/fL1fL1Zg/",
      plotsbasedir = "/work-zfs/lhc/heshy/LHEanomalouscouplings/fL1fL1Zg/plots/",
      lxplususername = "hroskes",
    )

for key, value in getconfiguration(socket.gethostname(), getpass.getuser()).iteritems():
  assert key not in globals()
  globals()[key] = value

repositorydir = os.path.realpath(repositorydir)

try:
    repositorydir
    plotsbasedir
    host
except NameError:
    raise ValueError("Who/where are you?\n{}\n{}".format(socket.gethostname(), getpass.getuser()))

LHE = True
usedata = False
showblinddistributions = False
unblinddistributions = False
unblindscans = False
useQGTagging = False
useVHMETTagged = True

usefastpdf = True
usefastpdfdouble = False

applym4lshapesystematicsUntagged = True
applym4lshapesystematicsVBFVHtagged = True
applym4lshapesystematicsggH = True
applym4lshapesystematicsggHUntagged = True
applym4lshapesystematicsdiagonal = True #VBF in VBFtagged, VH in VHHadrtagged
combinem4lshapesystematics = True

applyZXshapesystematicsVBFVHtagged = True
applyZXshapesystematicsUntagged = True
mergeZXVBFVHsystematics = True

applyJECshapesystematics = False
applyMINLOsystematics = True

assert unblindscans <= unblinddistributions <= showblinddistributions <= usedata <= (not LHE)

lumi2015 = 2.7

m4lmin, m4lmax = 105, 140

blindcut = lambda self: self.D_bkg() < 0.5

defaultnbins = 40



arrowsatminima = False
minimainlegend = True

if LHE:
    productionsforcombine = ["LHE_170509"]
    smearptelectron = 2.399/6
    smearptmuon = 2.169/6
    smearptjet = 18./6
else:
    productionsforcombine = ["170222"]

if len(productionsforcombine) == 1:
    productionforcombine = productionsforcombine[0]

