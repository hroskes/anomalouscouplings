import getpass
import os
import re
import socket


if (".cern.ch" in socket.gethostname() or "lxplus" in socket.gethostname()) and getpass.getuser() == "hroskes":
    host = "lxplus"
    repositorydir = "/afs/cern.ch/work/h/hroskes/anomalouscouplings_production/"
    plotsbasedir = "/afs/cern.ch/user/h/hroskes/www/anomalouscouplings_production/"
    svndir = "/afs/cern.ch/work/h/hroskes/AN/notes"

elif ("login-node" in socket.gethostname() or "compute" in socket.gethostname()) and getpass.getuser() == "jroskes1@jhu.edu":
    host = "MARCC"
    repositorydir = "/work-zfs/lhc/heshy/anomalouscouplings/"
    plotsbasedir = "/work-zfs/lhc/heshy/anomalouscouplings/plots/"

repositorydir = os.path.realpath(repositorydir)

try:
    repositorydir
    plotsbasedir
    host
except NameError:
    raise ValueError("Who/where are you?\n{}\n{}".format(socket.gethostname(), getpass.getuser()))

usedata = True
showblinddistributions = False
unblinddistributions = False
unblindscans = False
useQGTagging = False
useVHMETTagged = True

applym4lshapesystematics = True
combinem4lshapesystematics = True
applyZXshapesystematics = True
applyJECshapesystematics = False

assert unblindscans <= unblinddistributions <= showblinddistributions <= usedata

expectedscanluminosity = 35.867

m4lmin, m4lmax = 105, 140

blindcut = lambda self: self.D_bkg() < 0.5

productionsforcombine = ["170222"]

defaultnbins = 40
