import getpass
import os
import re
import socket

def getconfiguration(hostname, username):
  if (".cern.ch" in hostname or "lxplus" in hostname) and username == "hroskes":
    return dict(
      host = "lxplus",
      repositorydir = "/afs/cern.ch/work/h/hroskes/anomalouscouplings/",
      plotsbasedir = "/afs/cern.ch/user/h/hroskes/www/anomalouscouplings_multiparameter/" if hostname != socket.gethostname() else "/afs/cern.ch/user/h/hroskes/www/anomalouscouplings_HIG18002/",
      slcversion = 6,
      lxplususername = "hroskes",
      marccusername = "jroskes1@jhu.edu",
      svndir = "/afs/cern.ch/work/h/hroskes/AN/",
      connect = "hroskes@lxplus.cern.ch",
      name="heshy",
    )

  elif ("login-node" in hostname or "bc-login" in hostname or "compute" in hostname or "bigmem" in hostname) and username == "jroskes1@jhu.edu":
    return dict(
      host = "MARCC",
      repositorydir = "/work-zfs/lhc/heshy/anomalouscouplings/",
      repositorydir2015 = "/work-zfs/lhc/heshy/ICHEPanomalouscouplings/",
      plotsbasedir = "/work-zfs/lhc/heshy/anomalouscouplings/plots/",
      scratchdir = os.path.join("/scratch/users/", username, "tmparea", ""),
      slcversion = 7,
      lxplususername = "hroskes",
      marccusername = "jroskes1@jhu.edu",
      email = "heshyr@gmail.com",
      connect = "jroskes1@jhu.edu@gateway2.marcc.jhu.edu",
      name="heshy",
    )
  elif (".cern.ch" in hostname or "lxplus" in hostname) and username == "skyriaco":
    return dict(
      host = "lxplus",
      repositorydir = "/afs/cern.ch/work/s/skyriac2/anomalouscouplings/",
      plotsbasedir = "/afs/cern.ch/work/s/skyriac2/",
      slcversion = 6,
      lxplususername = "skyriaco",
      marccusername = "skyriac2@jhu.edu",
      svndir = "/afs/cern.ch/work/s/skyriaco/",
      connect = "skyriaco@lxplus.cern.ch",
      name="savvas",
    )

  elif ("login-node" in hostname or "bc-login" in hostname or "compute" in hostname or "bigmem" in hostname) and username == "skyriac2@jhu.edu":
    return dict(
      host = "MARCC",
      repositorydir = "/work-zfs/lhc/skyriaco/Combine/ggf/anomalouscouplings/",
      repositorydir2015 = "/work-zfs/lhc/skyriaco/Combine/ggf/anomalouscouplings/",
      plotsbasedir = "/work-zfs/lhc/skyriaco/Combine/ggf/anomalouscouplings/plots/",
      scratchdir = os.path.join("/scratch/users/", username, "tmparea", ""),
      slcversion = 7,
      lxplususername = "skyriaco",
      marccusername = "skyriac2@jhu.edu",
      email = "savasphy@gmail.com",
      connect = "skyriac2@jhu.edu@gateway2.marcc.jhu.edu",
      name="savvas",
    )

configuration = getconfiguration(socket.gethostname(), getpass.getuser())
for key, value in configuration.iteritems():
  assert key not in globals()
  globals()[key] = value

repositorydir = os.path.realpath(repositorydir)

try:
    repositorydir
    plotsbasedir
    host
except NameError:
    raise ValueError("Who/where are you?\n{}\n{}".format(socket.gethostname(), getpass.getuser()))

useQGTagging = False
useVHMETTagged = False

applym4lshapesystematics = True
applyZXshapesystematics = False
applyJECshapesystematics = True
applyMINLOsystematics = False
applySTXSsystematics = True

usebinbybin = False
staticmaxbins = 1000

lumi2015 = 2.7

m4lmin, m4lmax = 105, 140

blindcut = lambda self: self.D_bkg() < 0.5

defaultnbins = 40




arrowsatminima = False
minimainlegend = True

LHEsmearptelectron = 2.399/6
LHEsmearptmuon = 2.169/6
LHEsmearptjet = 18./6

usedata = True
showblinddistributions = True
unblinddistributions = True
unblindscans = True
productionsforcombine = ["190821_2016", "190821_2017", "190821_2018"]
#productionsforcombine = ["GEN_190908"]
separateZZWWVBFweights = True

assert unblindscans <= unblinddistributions <= showblinddistributions <= usedata

if len(productionsforcombine) == 1:
    productionforcombine = productionsforcombine[0]
