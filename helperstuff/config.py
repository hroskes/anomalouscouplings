import getpass
import os
import re
import socket

def getconfiguration(hostname, username):
  if (".cern.ch" in hostname or "lxplus" in hostname) and username == "hroskes":
    return dict(
      host = "lxplus",
      repositorydir = "/afs/cern.ch/work/h/hroskes/anomalouscouplings/",
      plotsbasedir = "/afs/cern.ch/user/h/hroskes/www/anomalouscouplings_HIG18002/",
      svndir = "/afs/cern.ch/work/h/hroskes/AN/",
      connect = "hroskes@lxplus.cern.ch",
    )

  elif ("login-node" in hostname or "compute" in hostname or "bigmem" in hostname) and username == "jroskes1@jhu.edu":
    return dict(
      host = "MARCC",
      repositorydir = "/work-zfs/lhc/heshy/anomalouscouplings/",
      repositorydir2015 = "/work-zfs/lhc/heshy/ICHEPanomalouscouplings/",
      plotsbasedir = "/work-zfs/lhc/heshy/anomalouscouplings/plots/",
      lxplususername = "hroskes",
      email = "heshyr@gmail.com",
      connect = "jroskes1@jhu.edu@gateway2.marcc.jhu.edu",
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

LHE = False
usedata = True
showblinddistributions = True
unblinddistributions = True
unblindscans = True

useQGTagging = False
useVHMETTagged = True

usefastpdf = True
usefastpdfdouble = False

applym4lshapesystematicsUntagged = True
applym4lshapesystematicsVBFVHtagged = True
applym4lshapesystematicsggH = False
applym4lshapesystematicsggHUntagged = False
applym4lshapesystematicsdiagonal = False #VBF in VBFtagged, VH in VHHadrtagged

getm4lsystsfromggHUntagged = False
getm4lsystsfromggH = True

combinem4lshapesystematics = False

applyZXshapesystematicsVBFVHtagged = True
applyZXshapesystematicsUntagged = True
mergeZXVBFVHsystematics = False
usenewZXsystematics = True

applyJECshapesystematics = True
applyMINLOsystematics = True

assert unblindscans <= unblinddistributions <= showblinddistributions <= usedata <= (not LHE)
assert getm4lsystsfromggHUntagged <= getm4lsystsfromggH

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
    productionsforcombine = ["180721_2016", "180721_2017"]

if len(productionsforcombine) == 1:
    productionforcombine = productionsforcombine[0]

